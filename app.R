#### AutobasisApp #####

###########################################################################
# LOAD REQUIRED LIBRARIES AND SET REPOSITORIES
###########################################################################

library(shiny)
library(shinyjs)
library(knitr)
library(dplyr)
library(R.utils)
library(data.table)
library(cupcake)

SNP.manifest <- copy(cupcake::SNP.manifest)
TITLE <- ' IMD Basis App '
## Header key
header_key <- c(chr = "CHR",
                pos = "POS",
                beta = "BETA",
                or = "OR",
                se = "SE",
                p_value = "P",
                p.value = "P",
                pval = "P",
                p = "P")

###########################################################################
# APP USER INTERFACE
###########################################################################

ui <- bootstrapPage(
  fluidPage( 
    column = 3, offset = 4, titlePanel(TITLE, windowTitle = TITLE)
    ),
  sidebarLayout(
    sidebarPanel(
        wellPanel(p("This shiny application allows you to project your GWAS summary statistics onto a lower dimensional basis that summarises the genetic architectures of 13 common immune-mediated diseases. For a more detailed explanation please see the 'help' tab on the right hand side.")),
        wellPanel(h3("Upload your data"),
        helpText("Upload files in .txt, .csv, .tsv. Compressed files (e.g. .txt.gz) are also accepted"),
        fileInput("user_data", "", accept = c(".txt", ".csv", ".tsv", ".txt.gz", ".tsv.gz", ".csv.gz"), multiple = F),
        textInput("trait_name", "Trait name", "IP-10"),
        textOutput("checkfile")),
        selectInput("PCDelta","Principal component", paste('PC',1:13,sep=''),selected='PC1',width="30%"),
        wellPanel(h3("Uploaded data overview"),tableOutput("QCtable"))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Forest Plot", wellPanel(verticalLayout(plotOutput("delta"),downloadButton(outputId = "downloadDelta",label = "Download Plot")))),
        tabPanel("Table", wellPanel(
                                    verticalLayout(tableOutput("table"),
                                        flowLayout(
                                        downloadButton(outputId = "downloadTable",label = "Download table"),
                                        checkboxInput(inputId = "downloadfulldata", label = "Include basis traits in download", value = FALSE))
                                    ))), 
        tabPanel("Help", column(12, includeMarkdown("help.Rmd")))
      )
    )
  )
)


###########################################################################
# SERVER SIDE
###########################################################################

server <- function(input, output, session) {
  
    ### Set file size limit. By default the file size limit is 5MB. Here limit is 200MB.
      options(shiny.maxRequestSize = 200*1024^2)
     ### Set reactive values
      rv <- reactiveValues(
        lastselectstage=NULL,
        selectstage='stageuploaddata'
      )
      PCDelta <- reactive({
        input$PCDelta
      })

###########################################################################
# HELPER FUNCTIONS
###########################################################################
      
    ### g.complement, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
      g.complement <- function (x) {
        x <- toupper(x)
        switches <- c(A = "t", T = "a", C = "g", G = "c")
        for (i in seq_along(switches)) x <- sub(names(switches)[i], switches[i], x)
        toupper(x)
      }
    
    ### g.rev, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
      g.rev <- function (x, sep = "/") {
        sapply(strsplit(x, sep), function(g) paste(rev(g), collapse = "/"))
      }
    
    ### g.class, as found in annotSnpStats package (github.com/chr1swallace/annotSnpStats/)
      g.class <- function (x, y) {
        if (!identical(names(x), names(y))) 
          stop("x and y must relate to same SNPs")
        mat <- matrix(FALSE, length(x), 4, dimnames = list(names(x), c("nochange", "rev", "comp", "revcomp")))
        mat[, "nochange"] <- x == y
        mat[, "rev"] <- x == g.rev(y)
        mat[, "comp"] <- x == g.complement(y)
        mat[, "revcomp"] <- x == g.rev(g.complement(y))
        indels <- x %in% c("I/D", "D/I")
        if (any(indels)) 
          mat[indels, c("comp", "revcomp")] <- FALSE
        ret <- character(nrow(mat))
        rs <- rowSums(mat)
        if (length(wh <- which(rs > 1))) 
          ret[wh] <- "ambig"
        if (length(wh <- which(rs == 0))) 
          ret[wh] <- "impossible"
        if (length(wh <- which(rs == 1))) 
          ret[wh] <- colnames(mat)[apply(mat[wh, , drop = FALSE], 1, which)]
        return(ret)
      }
    
    ### Check if there are missing columns
      missing_cols <- function(x){
      valid_headers <- c("CHR", "POS", "REF", "ALT", "SE", "P", "BETA")
      if(!all(valid_headers %in% names(x))){
        return(FALSE)
      }
      return(TRUE)
  
    }
    
    ### Check if data was correctly aligned to SNP manifest
     correct_alignment <- function(x){
      if(any(g.class(x$alleles.manifest, x$alleles) != "nochange")){
        "Some alleles do not match the SNP manifest file, this may be caused by data being in a different genome build, or by alleles needing some flipping. Please check."
      }else{
        NULL
      }
     }
     
###########################################################################
# Build the main data object
###########################################################################
     
     M <- reactive({
      if(is.null(input$user_data)){
        uploaded_data <- fread("data/Sample_dataset_B004_Ahola-Olli_27989323_1.tsv")
      } else {
        uploaded_data <- fread(input$user_data$datapath)
      }
      names(uploaded_data) <- header_key[names(uploaded_data)]
      # Check headers
      if(!"BETA" %in% names(uploaded_data) && "OR" %in% names(uploaded_data)){
        uploaded_data$BETA <- log(uploaded_data$OR)
      }
      shiny::validate(
        need(missing_cols(uploaded_data),"Required columns are missing in the uploaded file, please check the help tab for further details")
      )
      # Add two columns to facilitate alignment
      uploaded_data[,alleles:=paste(REF,ALT,sep="/")][,pid:=paste(CHR,POS,sep=":")]
      SNP.manifest[,alleles:=paste(ref_a1,ref_a2,sep="/")]
      
###########################################################################
# Align data with the manifest
###########################################################################
      
      M <- merge(uploaded_data,SNP.manifest[,.(pid,alleles)], by='pid', suffixes=c("",".manifest"))
      shiny::validate(
        correct_alignment(M)
      )
      return(M)
    })

###########################################################################
# Project data onto the basis
###########################################################################

    projected.userdata <- function(){
        projected.userdata <- cupcake::project_sparse(beta=M()$BETA, seb=M()$SE, pid=M()$pid)[,trait:=input$trait_name][]
        projected.userdata[, proj:=NULL] # Removed proj variable
        setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
        setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait")) # Removed Proj from here too
        return(projected.userdata)
    }
  
  # We extract the projected PCs and combine them with basis PC matrix, so we have a PC matrix that we can use for plotting
    combined.pcs <- function(){
        userdata.pcs <- projected.userdata()$Delta # Replaced proj by delta!
        combined.pcs <- rbind(userdata.pcs, basis$x[,1:13])
        rownames(combined.pcs)[1] <- input$trait_name
        return(combined.pcs)
    }

  # For Delta Plots, we'll calculate the significant PCs for the projected data
  # Then we'll need to project the basis traits onto the basis to characterize their PCs
  # Finally, we combine this data.table with the projected.userdata
   combined.deltaplot.dt <- function(){
       basistable.projected <- copy(cupcake::basis.trait.proj)
       setnames(basistable.projected, c("delta", "trait"), c("Delta", "Trait"))
       combined.deltaplot.dt <- rbind(basistable.projected[,.(PC,Delta,Var.Delta=0,Trait)],projected.userdata()[,.(PC,Delta,Var.Delta,Trait)])
       combined.deltaplot.dt[Var.Delta!=0, ci:=sqrt(Var.Delta) * 1.96]
       return(combined.deltaplot.dt)
   }
   

   # Prepare delta plot
   deltaplot.dt <- function(){
     deltaplot.dt <- combined.deltaplot.dt()[PC==PCDelta(),][order(Delta,decreasing = TRUE),]
     ndelta <- nrow(deltaplot.dt)
     ## dynamically compute x axis range
     idx <- which(!is.na(deltaplot.dt$ci))
     deltaplot.dt[is.na(ci),ci:=0]
     xl <- sort(c(deltaplot.dt$Delta[1] + deltaplot.dt$ci[1],deltaplot.dt$Delta[ndelta] - deltaplot.dt$ci[ndelta]) * 1.01)
     cols <- rep('black',nrow(deltaplot.dt))
     cols[idx] <- 'red'
     dotchart(deltaplot.dt$Delta,labels=deltaplot.dt$Trait,xlim=xl,pch=19, main=paste("Delta Plot", PCDelta()),xlab="Delta PC score",
              col=cols)
     ## add 95% confidence intervals
     with(deltaplot.dt[idx,],arrows(Delta-ci, idx, Delta+ci, idx, length=0.05, angle=90, code=3,col='red'))
     abline(v=0,col='red',lty=2)
   }

   # Prepare full table for download  
   full_data <- function(){
    full_data <- merge(projected.userdata()[, c("PC","P", "Trait")], combined.deltaplot.dt(), by = c("PC","Trait"), all.y = TRUE)
    setcolorder(full_data, c("PC", "Trait", "Delta", "Var.Delta", "P"))
    full_data[order(Trait,PC)]
    #setorder(full_data, "Trait")
   }
   
   
   
###########################################################################
# Graphic Output
###########################################################################


  ### Display PC table
  output$table <- renderTable({
      table_display <-  projected.userdata()[,.(PC,Var.Delta,Delta,P,Trait)]
      #table_display <- data.frame(PC = projected.userdata()$PC, Var.Delta = projected.userdata()$Var.Delta, Delta = projected.userdata()$Delta, P = projected.userdata()$P, Trait = projected.userdata()$Trait, stringsAsFactors = F) %>% arrange(P)
      table_display$Var.Delta <- sprintf("%.2e", table_display$Var.Delta)
      table_display$Delta <- sprintf("%.2e", table_display$Delta)
      table_display$P <- sprintf("%.2e", table_display$P)
      table_display
  }, width = '100%', align = 'r'
    )
  
  ### Download PC table
  output$downloadTable <- downloadHandler(
      filename = c('PCDelta_table.csv'),
      content = function(file) {
       if(input$downloadfulldata){
          write.csv(full_data(), file, row.names = FALSE, quote = FALSE)
        }else{
          write.csv(projected.userdata()[, c("PC", "Var.Delta", "Delta", "P", "Trait")], file, row.names=FALSE, quote = FALSE)
        }
     }
  )
  ### Display QC table
  output$QCtable <- renderTable({
    vars <- c("SNPs in basis","Uploaded SNPs overlapping basis", "Overall P-value", "Most significant component")
    n.snps.basis <- nrow(SNP.manifest)
    n.snps <- nrow(M())
    SNP_overlap <- sprintf("%d (%.1f%%)",n.snps,(n.snps/n.snps.basis)*100)
    overall_p <- sprintf("%1.0e",projected.userdata()$p.overall[1])
    minc.idx <- which.min(projected.userdata()$P)
    min.component_txt <- sprintf("%s (%1.0e)",projected.userdata()$PC[minc.idx],projected.userdata()$P[minc.idx])
    values <- c(n.snps.basis,SNP_overlap, overall_p, min.component_txt)
    QCTable <- data.frame(vars, values)
    QCTable
  }, colnames = FALSE)
  
  
  ### Display Delta plot
  output$delta <- renderPlot({
      deltaplot.dt()
  },execOnResize = TRUE)
  
  ### Download Delta plot 
  output$downloadDelta <- downloadHandler(
    filename = function(){
      paste(input$trait_name, input$PCDelta, "deltaplot.pdf", sep = "_")
    },
    content = function(file){
      pdf(file)
      deltaplot.dt()
      dev.off()
    }
  )
    
}

shinyApp(ui, server)
