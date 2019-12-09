#### AutobasisApp #####

###########################################################################
# LOAD REQUIRED LIBRARIES AND SET REPOSITORIES
###########################################################################

library(shiny)
library(shinyjs)
library(BiocManager)
options(repos = BiocManager::repositories())
library(dplyr)
library(knitr)
library(R.utils)
library(data.table)
library(cupcake)

###########################################################################
# LOAD DATASETS
###########################################################################

data(burren_imd_13_traits)

SPARSE_BASIS_EXTDATA <- system.file("extdata/sparse_imd_basis", package="cupcake")
SNP_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_shrinkage.RDS')
BASIS_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_basis.RDS')
GWAS_DATA_DIR <- file.path(SPARSE_BASIS_EXTDATA,'/gwas_data/')
basis.gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR)
basis <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)

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

ui <- fluidPage(

  headerPanel(title = "", windowTitle = "AutobasisApp"),
  wellPanel(
    titlePanel(
      fluidRow(
        column(10,
               div("AutobasisApp", style = "margin-left:10px;"),
               div("Shiny Web Application for projecting user GWAS data onto a 13 Immune-mediated traits Basis", style = "padding:10px; font-size:70%;"))
              ))),

  
###########################################################################
# STAGEHELP
###########################################################################
  
  conditionalPanel(condition = "input.selectstage == 'stagehelp'",
                   
                   div(style="display: inline-block;vertical-align:middle;", h3("Help")),
                   div(style="display: inline-block;vertical-align:middle;", h3(" ")),
                   div(style="display: inline-block;vertical-align:middle;", actionButton("returnButton", icon = icon("reply"),label="")),
                   fluidRow(column(12,HTML("<br>"))),
                   
                   tabsetPanel(
                     tabPanel("About", column(12, includeMarkdown("about.md"))),
                     
                     textInput("selectstage", label="", value="stageuploaddata", width='300px') # For some reason this line is crucial for the whole app to work properly
                   )),
      h3(" "),

###########################################################################
# STAGE OUTPUT
###########################################################################
  
        conditionalPanel(condition = "input.selectstage == 'stageuploaddata'",
            div(style="display: inline-block;vertical-align:middle;", actionButton("helpButton", strong("Help") )),
            h2(""),
            sidebarLayout(                 
            sidebarPanel(
              h3("Upload data"),
              helpText("Upload files in .txt, .csv, .tsv. Compressed files (e.g. .txt.gz) are also accepted"),
              fileInput("user_data", "", accept = c(".txt", ".csv", ".tsv", ".txt.gz", ".tsv.gz", ".csv.gz"), multiple = F),
              helpText("Provide the name of the trait you want to project onto the basis"),
              textInput("trait_name", "Trait name", "IP-10"),
              sliderInput("PCDelta", "Principal component (Delta Plot)", 1, 13, 1, pre = "PC")
            ),
            mainPanel(
              fluidRow(
                column(9, align = "center", offset = 1,
                       plotOutput("delta"))
              ),
              fluidRow(
                downloadButton(outputId = "downloadDelta",label = "Download Delta Plot")
              ),
              fluidRow(
                h4("Quality Control"),
                tableOutput("QCtable")
              ),
              fluidRow(
                h4("Projected data results"),
                tableOutput("table")
              ),
              fluidRow(
                  checkboxInput(inputId = "downloadfulldata", label = "I want full results", value = FALSE)
              ),
              fluidRow(
                downloadButton(outputId = "downloadTable",label = "Download table")
                
              )
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
  
    ### Switch between different pannels
      observeEvent(input$helpButton , {
        rv$lastselectstage<-input$selectstage
        updateTextInput(session, "selectstage", value = 'stagehelp')
        toggle("helpButton")
      })
      
      observeEvent(input$returnButton , {
        updateTextInput(session, "selectstage", value = rv$lastselectstage)
        toggle("helpButton")
      })
  
     ### Set reactive values
      rv <- reactiveValues(
        lastselectstage=NULL,
        selectstage='stageuploaddata'
      )
  
    ### Transform initial slider values to reactives
      PCScatterplot <- reactive({
        paste("PC", input$PC, sep = "")
      })
      
      PCDelta <- reactive({
        paste("PC", input$PCDelta, sep = "")
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
        missing_cols <- valid_headers[!valid_headers %in% names(x)]
        print("These column names seem to be missing: ", cat(valid_headers[!valid_headers %in% names(x)], sep = ", "), ". Please provide valid headers.")
      }
      else{
        NULL
      }
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
        missing_cols(uploaded_data)
      )
      # Add two columns to facilitate alignment
      uploaded_data[,alleles:=paste(REF,ALT,sep="/")][,pid:=paste(CHR,POS,sep=":")]
      SNP.manifest[,alleles:=paste(ref_a1,ref_a2,sep="/")]
      
###########################################################################
# Align data with the manifest
###########################################################################
      
      M <- merge(uploaded_data,SNP.manifest[,.(pid,alleles)], by='pid', suffixes=c("",".manifest"))
      # SNP missing in user data
      # if(nrow(SNP.manifest) > nrow(M)) warning("There are ", nrow(SNP.manifest) - nrow(M), " less SNPs in data than in base. Please check.")
      # # NOT SURE ABOUT WHAT TO DO WHEN THIS DOES NOT EQUAL ZERO
      
      shiny::validate(
        correct_alignment(M)
      )
      return(M)
    })

###########################################################################
# Project data onto the basis
###########################################################################

    projected.userdata <- reactive({
        projected.userdata <- cupcake::project_sparse(beta=M()$BETA, seb=M()$SE, pid=M()$pid)[,trait:=input$trait_name][]
        projected.userdata[, proj:=NULL] # Removed proj variable
        setnames(projected.userdata, c("var.proj", "delta", "p", "trait"), c("Var.Delta", "Delta", "P", "Trait")) # Changed var.proj to var.delta
        setcolorder(projected.userdata, c("PC", "Var.Delta", "Delta", "p.overall", "z", "P", "Trait")) # Removed Proj from here too
        return(projected.userdata)
    })
  
  # We extract the projected PCs and combine them with basis PC matrix, so we have a PC matrix that we can use for plotting
    combined.pcs <- reactive({
        userdata.pcs <- projected.userdata()$Delta # Replaced proj by delta!
        combined.pcs <- rbind(userdata.pcs, basis$x[,1:13])
        rownames(combined.pcs)[1] <- input$trait_name
        return(combined.pcs)
    })

  # For Delta Plots, we'll calculate the significant PCs for the projected data
  # Then we'll need to project the basis traits onto the basis to characterize their PCs
  # Finally, we combine this data.table with the projected.userdata
   combined.deltaplot.dt <- reactive({
       basis.gwas.DT[,c('beta','seb'):=list(log(or),1) ]
       basistable.projected <- lapply(split(basis.gwas.DT,basis.gwas.DT$trait),function(x){
         tt <- x$trait %>% unique
         cupcake::project_sparse(beta=x$beta,seb=x$seb,pids=x$pid)[,trait:=tt]
       }) %>% rbindlist
       setnames(basistable.projected, c("var.proj", "delta", "trait"), c("Var.Delta", "Delta", "Trait"))
       combined.deltaplot.dt <- rbind(basistable.projected[,.(PC,Delta,Var.Delta=0,Trait)],projected.userdata()[,.(PC,Delta,Var.Delta,Trait)])
       combined.deltaplot.dt[Var.Delta!=0, ci:=sqrt(Var.Delta) * 1.96]
       return(combined.deltaplot.dt)
   })
   

   # Prepare delta plot
   deltaplot.dt <- reactive({
     deltaplot.dt <- combined.deltaplot.dt()[PC==PCDelta(),][order(Delta,decreasing = TRUE),]
     idx <- which(!is.na(deltaplot.dt$ci))
     cols <- rep('black',nrow(deltaplot.dt))
     cols[idx] <- 'red'
     dotchart(deltaplot.dt$Delta,labels=deltaplot.dt$Trait,xlim=c(-0.1,0.05),pch=19, main=paste("Delta Plot", PCDelta()),xlab="Delta PC score",
              col=cols)
     ## add 95% confidence intervals
     with(deltaplot.dt[idx,],arrows(Delta-ci, idx, Delta+ci, idx, length=0.05, angle=90, code=3,col='red'))
     abline(v=0,col='red',lty=2)
   })

   # Prepare full table for download  
   ## Full dataset
   full_data <- reactive({
    full_data <- merge(projected.userdata()[, c("PC","P", "Trait")], combined.deltaplot.dt(), by = c("PC","Trait"), all.y = TRUE)
    setcolorder(full_data, c("PC", "Trait", "Delta", "Var.Delta", "P"))
    setorder(full_data, "Trait")
   })
   
   
   
###########################################################################
# Graphic Output
###########################################################################


  ### Display PC table
  output$table <- renderTable({
      table_display <- data.frame(PC = projected.userdata()$PC, Var.Delta = projected.userdata()$Var.Delta, Delta = projected.userdata()$Delta, P = projected.userdata()$P, Trait = projected.userdata()$Trait, stringsAsFactors = F) %>% arrange(P)
      table_display$Var.Delta <- sprintf("%.4e", table_display$Var.Delta)
      table_display$Delta <- sprintf("%.4e", table_display$Delta)
      table_display$P <- sprintf("%.2e", table_display$P)
      table_display
  }, width = '100%'
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
    vars <- c("Data v. Manifest overlapping SNPs", "Data v. Manifest overlapping SNPs (%)", "Overall P-value")
    SNP_overlap <- sprintf("%d",nrow(M()))
    SNP_overlap_percentage <- sprintf("%.2f",(nrow(M())/nrow(SNP.manifest))*100)
    overall_p <- sprintf("%7.2e",projected.userdata()$p.overall[1])
    values <- c(SNP_overlap, SNP_overlap_percentage, overall_p)
    QCTable <- data.frame(vars, values)
    QCTable
  }, width = "50%", colnames = FALSE)
  
  # ### Display PCA Scatterplot
  # output$PCAScatterplot <- renderPlot({
  #     scatterplot.dt <- data.table(traits = rownames(combined.pcs()), combined.pcs())
  #     xlim <- c(min(scatterplot.dt[, get(PCScatterplot()[1])]),max(scatterplot.dt[, get(PCScatterplot()[1])])) * 1.2
  #     ylim <- c(min(scatterplot.dt[, get(PCScatterplot()[2])]),max(scatterplot.dt[, get(PCScatterplot()[2])])) * 1.2
  #     cols <- rep('black', nrow(scatterplot.dt))
  #     cols[scatterplot.dt$traits == input$trait_name] <- 'red'
  #     plot(scatterplot.dt[, get(PCScatterplot()[1])],scatterplot.dt[, get(PCScatterplot()[2])],type='n',xlim=xlim,ylim=ylim,main="PCA scatterplot", 
  #          xlab = PCScatterplot()[1], ylab = PCScatterplot()[2], col=cols)
  #     abline(h=0, v=0, col="red", lty = 2, lwd = 1)
  #     text(scatterplot.dt[, get(PCScatterplot()[1])], scatterplot.dt[, get(PCScatterplot()[2])],labels=scatterplot.dt$traits, cex= 0.8, adj=c(0,0), col=cols)
  #     points(scatterplot.dt[, get(PCScatterplot()[1])], scatterplot.dt[, get(PCScatterplot()[2])],cex=0.5,pch=19, col=cols)
  # })
  
  ### Display Delta plot
  output$delta <- renderPlot({
      deltaplot.dt()
  })
  
  ### Download Delta plot 
  output$downloadDelta <- downloadHandler(
    filename = function(){
      paste(input$trait_name, input$PCDelta, "deltaplot.pdf", sep = "_")
      },
    content = function(file){
      pdf(file)
      # Had to draw the plot twice, since using the reactive plot didn't work
      deltaplot.dt <- combined.deltaplot.dt()[PC==PCDelta(),][order(delta,decreasing = TRUE),]
      idx <- which(!is.na(deltaplot.dt$ci))
      cols <- rep('black',nrow(deltaplot.dt))
      cols[idx] <- 'red'
      dotchart(deltaplot.dt$delta,labels=deltaplot.dt$trait,xlim=c(-0.1,0.05),pch=19, main=paste("Delta Plot", PCDelta()),xlab="Delta PC score",
               col=cols)
      ## add 95% confidence intervals
      with(deltaplot.dt[idx,],arrows(delta-ci, idx, delta+ci, idx, length=0.05, angle=90, code=3,col='red'))
      abline(v=0,col='red',lty=2)
      dev.off()
    }
  )
    
}

shinyApp(ui, server)
