#### Shiny-basis project #####

#### LOAD REQUIRED LIBRARIES AND SET REPOSITORIES ####
library(shiny)
library(BiocManager)
options(repos = BiocManager::repositories())
library(dplyr)
library(R.utils)
library(data.table)
library(DT)
library(cupcake)
library(annotSnpStats) # on github.com/chr1swallace/annotSnpStats
# library(cowplot)
library(RColorBrewer)
# library(ggplot2)

#### LOAD DATASETS ####
data(burren_imd_13_traits)

## Dummy input
# input <- list(user_data = fread("EOC_Kanai_29403010.txt.gz"), trait_name = "EOC", PC=c(1,2), PCDelta = 1)

SPARSE_BASIS_EXTDATA <- system.file("extdata/sparse_imd_basis", package="cupcake")
SNP_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_snp_manifest.RDS')
TRAIT_MANIFEST_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/sparse_trait_manifest.tab')
SHRINKAGE_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_shrinkage.RDS')
BASIS_FILE <- file.path(SPARSE_BASIS_EXTDATA,'support/13_trait_sparse_basis.RDS')
GWAS_DATA_DIR <- file.path(SPARSE_BASIS_EXTDATA,'/gwas_data/')
basis.gwas.DT<-get_gwas_data(TRAIT_MANIFEST_FILE,SNP_MANIFEST_FILE,GWAS_DATA_DIR)
basis <- readRDS(BASIS_FILE)
shrink.DT <- readRDS(SHRINKAGE_FILE)
# Header key
header_key <- c(chr = "CHR",
                pos = "POS",
                beta = "BETA",
                or = "OR",
                se = "SE",
                p_value = "P",
                pval = "P",
                p = "P")



#### APP USER INTERFACE ####
ui <- fluidPage(

  headerPanel(title = "", windowTitle = "AutoBasisApp"),
  wellPanel(
    titlePanel(
      fluidRow(
        column(10,
               div("AutoBasisApp", style = "margin-left:10px;"),
               div("Shiny Web Application for projecting user GWAS data onto a 13 Immune-mediated traits Basis", style = "padding:10px; font-size:70%;"))
              ))),
  sidebarLayout(
        sidebarPanel(
          h3("Upload data"),
          helpText("Upload files in .txt, .csv, .tsv. Compressed files (e.g. .txt.gz) are also accepted"),
          fileInput("user_data", "", accept = c(".txt", ".csv", ".tsv", ".txt.gz", ".tsv.gz", ".csv.gz"), multiple = F),
          helpText("Provide the name of the trait you want to project onto the basis"),
          textInput("trait_name", "Trait name", "User trait"),
          sliderInput("PC", "Principal components (Scatter Plot)", 1, 13, c(1, 2), pre = "PC"),
          sliderInput("PCDelta", "Principal component (Delta Plot)", 1, 13, 1, pre = "PC"),
          actionButton("startAnalysisButton",strong( "Start analysis"), style='background-color: #F1A8A8')
        ),
        mainPanel(
          fluidRow(
            splitLayout(cellWidths = c("50%", "50%"), plotOutput("PCAScatterplot"), plotOutput("forest"))
          ),
          fluidRow(
              splitLayout(cellWidths = c("50%", "50%"), DTOutput("table"), plotOutput("delta"))
          )
          
        )
        )
)

#### SERVER SIDE CODE ####
server <- function(input, output, session) {
  
  ## (1) SANITY CHECKS ON INPUT FILE
  # By default the file size limit is 5MB. Here limit is 70MB.
   options(shiny.maxRequestSize = 200*1024^2)
  
  # Increase memory limit
  # memory.size(max = FALSE)
  observeEvent(input$startAnalysisButton, {
    
    PCScatterplot <- reactive({
      paste("PC", input$PC, sep = "")
    })
    
    PCDelta <- reactive({
      paste("PC", input$PCDelta, sep = "")
    })
    
    M <- reactive({
        uploaded_data <- fread(input$user_data$datapath)
        names(uploaded_data) <- header_key[names(uploaded_data)]
        # Check headers
        if(!"BETA" %in% names(uploaded_data) && "OR" %in% names(uploaded_data)){
           uploaded_data$BETA <- log(uploaded_data$OR)
        }
        valid_headers <- c("CHR", "POS", "REF", "ALT", "SE", "P", "BETA")
        if(!all(valid_headers %in% names(uploaded_data))){
          missing_cols <- valid_headers[!valid_headers %in% names(uploaded_data)]
          stop("These column names seem to be missing: ", cat(valid_headers[!valid_headers %in% names(uploaded_data)], sep = ", "), ". Please provide valid headers.")
        }
        
        # Add two columns to facilitate alignment
        uploaded_data[,alleles:=paste(REF,ALT,sep="/")][,pid:=paste(CHR,POS,sep=":")]
        SNP.manifest[,alleles:=paste(ref_a1,ref_a2,sep="/")]
        
        ## (2) ALIGN THE DATA TO SNP MANIFEST
        
        M <- merge(uploaded_data,SNP.manifest[,.(pid,alleles)], by='pid', suffixes=c("",".manifest"))
        # SNP missing in user data
        if(nrow(SNP.manifest) > nrow(M)) warning("There are ", nrow(SNP.manifest) - nrow(M), " less SNPs in data than in base. Please check.")
        # NOT SURE ABOUT WHAT TO DO WHEN THIS DOES NOT EQUAL ZERO
        
        # Sanity check for alignment
        if(any(g.class(M$alleles.manifest, M$alleles) != "nochange")){
          ### EXTRA HARMONIZING STEPS. TO BE DEVELOPED
          stop("Some alleles do not match the SNP manifest file, this may be caused by data being in a different genome build, or by alleles needing some flipping. Please check.")
        }
        return(M)
    })


  ## (3) PROJECT ONTO THE BASIS

    projected.userdata <- reactive({
        projected.userdata <- cupcake::project_sparse(beta=M()$BETA, seb=M()$SE, pid=M()$pid)[,trait:=input$trait_name][]
        projected.userdata[, p_adj := p.adjust(p, method = "bonferroni"),]
        setcolorder(projected.userdata, c("PC", "proj", "var.proj", "delta", "p.overall", "z", "p", "p_adj", "trait"))
        return(projected.userdata)
    })
  
  # We extract the projected PCs and combine them with basis PC matrix, so we have a PC matrix that we can use for plotting
    combined.pcs <- reactive({
        userdata.pcs <- projected.userdata()$proj
        combined.pcs <- rbind(userdata.pcs, basis$x[,1:13])
        rownames(combined.pcs)[1] <- input$trait_name
        return(combined.pcs)
    })

  # For Delta Plots, we'll calculate the significant PCs for the projected data
  # Then we'll need to project the basis traits onto the basis to characterize their PCs
  # Finally, we combine this data.table with the projected.userdata
  basis.gwas.DT[,c('beta','seb'):=list(log(or),1) ]
  basistable.projected <- lapply(split(basis.gwas.DT,basis.gwas.DT$trait),function(x){
    tt <- x$trait %>% unique
    cupcake::project_sparse(beta=x$beta,seb=x$seb,pids=x$pid)[,trait:=tt]
  }) %>% rbindlist
  combined.deltaplot.dt <- rbind(basistable.projected[,.(PC,delta,var.proj=0,trait)],projected.userdata()[,.(PC,delta,var.proj,trait)])
  combined.deltaplot.dt[var.proj!=0,ci:=sqrt(var.proj) * 1.96]

  ## (4) GRAPHIC OUTPUT
  
      ## PC table
      output$table <- DT::renderDataTable({
      projected.userdata()
      })
    
      ## PCA Scatterplot
      output$PCAScatterplot <- renderPlot({
      scatterplot.dt <- data.table(traits = rownames(combined.pcs()), combined.pcs())
      xlim <- c(min(scatterplot.dt[, get(PCScatterplot()[1])]),max(scatterplot.dt[, get(PCScatterplot()[1])])) * 1.2
      ylim <- c(min(scatterplot.dt[, get(PCScatterplot()[2])]),max(scatterplot.dt[, get(PCScatterplot()[2])])) * 1.2
      cols <- rep('black', nrow(scatterplot.dt))
      cols[scatterplot.dt$traits == input$trait_name] <- 'red'
      plot(scatterplot.dt[, get(PCScatterplot()[1])],scatterplot.dt[, get(PCScatterplot()[2])],type='n',xlim=xlim,ylim=ylim,main="PCA scatterplot", 
           xlab = PCScatterplot()[1], ylab = PCScatterplot()[2], col=cols)
      text(scatterplot.dt[, get(PCScatterplot()[1])],scatterplot.dt[, get(PCScatterplot()[2])],labels=scatterplot.dt$traits, cex= 0.8, adj=c(0,0), col=cols)
      points(scatterplot.dt[, get(PCScatterplot()[1])],scatterplot.dt[, get(PCScatterplot()[2])],cex=0.5,pch=19, col=cols)
      })
      
      ## Delta plot
      
      output$delta <- renderPlot({
        deltaplot.dt <- combined.deltaplot.dt[PC==PCDelta(),][order(delta,decreasing = TRUE),]
        idx <- which(!is.na(deltaplot.dt$ci))
        cols <- rep('black',nrow(deltaplot.dt))
        cols[idx] <- 'red'
        dotchart(deltaplot.dt$delta,labels=deltaplot.dt$trait,xlim=c(-0.1,0.05),pch=19, main=paste("Delta Plot", PCDelta()),xlab="Delta PC score",
                                       col=cols)
          ## add 95% confidence intervals
        with(deltaplot.dt[idx,],arrows(delta-ci, idx, delta+ci, idx, length=0.05, angle=90, code=3,col='red'))
        abline(v=0,col='red',lty=2)
      })
      ## Forest plot
      output$forest <- renderPlot({
      combined.pcs() %>% dist %>% hclust %>% plot(., main = 'Forest Plot')
      })
  }) #ObserveEvent(button)
}




shinyApp(ui, server)