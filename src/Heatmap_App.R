#Load libraries

library("shiny") 
library("shinyjs")
library("d3heatmap") 
library("vecsets")
library("stringr")

##' Program to create a Heatmap from the following three input files with option to save the generated file and plot a d3heatmap
##' Input files 
##'          1. autism_prediction.v2.rscores.exp100.txt contains all the genes and their ASD-association scores.
##'          2. brainspan.spatio-temporal.rma.st-rowmodz-stoufferzge2.gmt contains the genes specific each spatiotemporal window.
##'          3. reg-per.sub.txt contains the subset of brain regions and developmental periods that need to be considered.
##' @author Sreeja Ashok

if (interactive()) {
  runApp(shinyApp(
    ui = fluidPage(
      shinyjs::useShinyjs(),
      sidebarLayout(
        sidebarPanel( 
          radioButtons("rbSelection", "Choose your option:",
                       choiceNames = list("Create a Heatmap" ,
                                          "Plot Heatmap"),
                       choiceValues = list("construct", "load")
          ),
          conditionalPanel(condition = "input.rbSelection == 'construct'",
                           fileInput( 'fileLoaded', 'Choose File',
                                      accept = c( 'text/csv', 
                                                  'text/comma-separated-values,text/plain', 
                                                  '.csv')),
                           selectInput("inSelect", "Select File", choices = ""),
                           checkboxInput("header", "Header", TRUE),
                           radioButtons(inputId = 'sep', label = 'Separator', choices = c( Comma=',',Semicolon=';',Tab='\t', 
                                                                                           Space=''), selected = ','),
                           actionButton("heatmapmatrix", "Build Heatmap Matrix"),
                           downloadButton("saveFile", "Save  Heatmap Matrix"),
                           div(style="display:inline-block;width:45%;text-align: center;", selectInput("palette", "Palette", 
                                                                                                       c("YlOrRd","RdYlBu", "Greens", "Blues"))),
                           div(style="display:inline-block;width:45%;text-align: center;",selectInput("dendrogram",
                                                                                                      "Dendogram",c("both", "row","column", "none"))),
                           shinyjs::hidden(p(id = "text1", "Processing..."))
          ),
          conditionalPanel(condition = "input.rbSelection == 'load'",
                           fileInput( 'heatMapfileLoaded', 'Choose File',
                                      accept = c( 'text/csv', 
                                                  'text/comma-separated-values,text/plain', 
                                                  '.csv')),
                           div(style="display:inline-block;width:45%;text-align: center;", selectInput("palette1", "Palette", 
                                                                                                       c("YlOrRd","RdYlBu", "Greens", "Blues"))),
                           div(style="display:inline-block;width:45%;text-align: center;", selectInput("dendrogram1", 
                                                                                                       "Dendogram",c("both", "row","column", "none")))
          ) # conditional panel 
        ), #sidebarpanel
        mainPanel(
          conditionalPanel(condition = "input.rbSelection == 'construct'",
                           tabsetPanel(
                             tabPanel("Data", dataTableOutput("contents")),
                             tabPanel("Summary", verbatimTextOutput("sum")),
                             tabPanel("Heatmap",fluidRow(d3heatmapOutput("heatmap",width = "100%", 
                                                                         height="500px")))
                           )
          ),
          conditionalPanel(condition = "input.rbSelection == 'load'",
                           fluidRow(d3heatmapOutput("heatmap1",width = "100%",height = "500px"))
                           #  fluidRow(plotlyOutput("heatmap1"))
          )
        ) # main panel
      ) # sidebarLayout
    ), #ui
    
    server = function(input, output, session) {
      
      # Variable to store the filenames and filepaths loaded
      fileNameList<<-NULL
      filePathList<<-NULL
      fileAutismData<<-NULL
      fileBrainData <<-NULL
      fileRegData <<-NULL
      HeatMapMatrix<<-NULL
      HeatMapMatrix1<<-NULL
      
      # Disable UI components initially
      shinyjs::disable("saveFile")
      shinyjs::disable("heatmapmatrix")
      shinyjs::disable("palette")
      shinyjs::disable("dendrogram")
      
      # Variable to set the view option of columnn data 
      dimSize <- 30
      # Display the contents of each file based on file selection
      output$contents <- renderDataTable({
        if( is.null(fileData())){return()}
        if(ncol(fileData()) > dimSize ){
          fileData()[1:dimSize]
        } else
          fileData()
      }) 
      
      # Build the heatmap matrix
      observeEvent(input$heatmapmatrix, {
        if ( (!is.null(fileAutismData)) && (!is.null(fileBrainData)) && (!is.null(fileRegData)) ) {
          buildHeatMapMatrix(fileAutismData, fileBrainData, fileRegData)
        } else {
          showNotification("Please load all the suitable files")
        }
      })
      
      # Save the Heatmap matrix generated. Mention the filename and extension as "filename.csv" in the dialog window 
      output$saveFile <- downloadHandler(
        filename = function() {	
          paste("heatmapMatrix",".csv",sep="")
        },
        content = function(file) {
          write.csv(HeatMapMatrix, file)
          html("text1", paste("Heatmap Matrix is saved"))
        }
      )
      
      # Function to disable the UI components
      disableUI <- function(){
        shinyjs::disable("heatmapmatrix")
        shinyjs::disable("rbSelection")
        shinyjs::disable("fileLoaded")
        shinyjs::disable("inSelect")
        shinyjs::disable("header")
        shinyjs::disable("saveFile")
        shinyjs::disable("Summary")
        shinyjs::disable("Data")
      }
      
      # Function to enable the UI components
      enableUI <- function(){
        shinyjs::enable("rbSelection")
        shinyjs::enable("fileLoaded")
        shinyjs::enable("inSelect")
        shinyjs::enable("header")
        shinyjs::enable("saveFile")
        shinyjs::enable("Summary")
        shinyjs::enable("Data")
        shinyjs::enable("palette")
        shinyjs::enable("dendrogram")
      }
      
      # Function to build the Heatmap matrix from the 3 files given
      buildHeatMapMatrix <- function(file1,file2,file3) {
        disableUI()
        shinyjs::show("text1")
        showNotification("Processing .......")
        fileAutismScores <- read.csv(file1,header=TRUE)
        fileBrainSpanData <- read.csv(file2,header=FALSE)
        fileReg <- read.csv(file3,header=TRUE)
        
        fileAutismScores <- as.matrix(fileAutismScores)
        fileBrainSpanData <- as.matrix(fileBrainSpanData)
        fileReg <- as.matrix(fileReg)
        geneSet_fileBrainSpanData <- fileBrainSpanData[,-1:-2]
        
        # Get all non zero values in fileBrainSpanData and convert it to vector 
        brain_Data <- as.vector(as.double(geneSet_fileBrainSpanData[!is.na(geneSet_fileBrainSpanData)]))
        aut_Scores <- as.vector(fileAutismScores[,1])
        in_aut_not_brain <- vsetdiff(aut_Scores,brain_Data)
        #in_brain_not_aut <- vsetdiff(brain_Data,aut_Scores) 
        autScores_Subset <- fileAutismScores [aut_Scores %in% brain_Data,  ]
        
        # Define the rows and colums of Heapmap matrix 
        
        rowHeatmap <- nrow(fileReg[fileReg[,3] == "reg",])
        colHeatmap  <-   nrow(fileReg[fileReg[,3] == "per",])
        HeatMapMatrix <- matrix( 0 ,nrow = rowHeatmap,ncol = colHeatmap)
        rownames(HeatMapMatrix)<- fileReg[fileReg[,3] == "reg",1]
        colnames(HeatMapMatrix)<- fileReg[fileReg[,3] == "per",2]
        
        # Initialize a one dimensional array(vector) of size 10,000 initialized with zero values
        N <- 10000
        aut_Score_randVect <- as.vector(rep(0, N))
        
        # Function to get the row index of HeatmapMatrix 
        setRowCounter <- function(x){
          i <- 0
          function(){
            i <<- i+1
            if(i> nrow(HeatMapMatrix))
              i <<- 1
            i   
          }
        }
        rowCounter <- setRowCounter()
        
        # Function to get the column index of HeatmapMatrix 
        setColCounter <- function(x){
          i <- 0
          j <- 0
          function(){
            i <<- i+1
            if((i-1) %% nrow(HeatMapMatrix) == 0)
              j <<- j+1
            j   
          }
        }
        colCounter <- setColCounter()
        
        # Computing the Q-val of the genes mapped in each cell based on information from brainspan-genesets and autism scores.
        HeatMapMatrix <-
          apply(HeatMapMatrix, c(1,2), function(x) {
            # Get the rowname+Colname of fileReg to access the genes associated with that cell from the file , BrainSpanData 
            location <-
              paste(rownames(HeatMapMatrix)[rowCounter()],str_pad(as.character(fileReg[fileReg[,2] ==
                                                                                         colnames(HeatMapMatrix)[colCounter()],1]),2,pad = "0"), sep = ".")
            # Extract the number of genes associated with the region
            geneCount <-
              strsplit (as.character(fileBrainSpanData [fileBrainSpanData[,1] == location,2]),'[()]')[[1]][2]
            # Extract all the genesIds mapped to a cell
            geneIDs <-
              fileBrainSpanData[which(fileBrainSpanData[,1] == location), 3:(as.numeric(geneCount) + 2)]
            # Compute the mean score
            mean_Score <-
              mean(autScores_Subset[autScores_Subset[,1] %in% as.integer(geneIDs),2])
            geneCount_AutScore <-
              as.numeric(nrow(autScores_Subset[autScores_Subset[,1] %in% as.integer(geneIDs),]))
            aut_Score_randVect <-
              sapply(aut_Score_randVect, function(x)
                mean(autScores_Subset[sample(nrow(autScores_Subset),
                                             size = geneCount_AutScore),2]))
            aut_Score_randVect[1] <- mean_Score
            p_Val <-
              (sum(aut_Score_randVect >= aut_Score_randVect[1]) / N)
            p_Adj <-  p.adjust(p_Val, "BH")
            - log10(p_Adj)
          })
        HeatMapMatrix <<- HeatMapMatrix
        enableUI()
        html("text1", paste("Heatmap Matrix is generated"))
      }
      
      # Display the summary of the file information
      output$sum <- renderPrint({
        if( is.null(fileData())){
          return("File not loaded")
        }
        summary(fileData())
      })
      
      # Plot  heatmap for  option 1
      output$heatmap <- renderD3heatmap({
        if (is.null(HeatMapMatrix)){
          print("HeatMapMatrix is null")
          return(NULL)
        }
        d3heatmap(
          HeatMapMatrix,
          colors = input$palette,
          xaxis_height = 200,
          yaxis_width = 120,
          digits = 4L,
          dendrogram = input$dendrogram 
        )
      })
      
      # Read the input file
      observeEvent(input$inSelect,{
        fileData()
      })
      
      # Read the input file
      fileData<- reactive({
        inFile <- filePathList[match(input$inSelect, fileNameList)]
        if ( is.null(inFile))
          return(NULL)
        read.csv(inFile, sep=input$sep, header=input$header)
      })
      
    # view heatmap for option2
       observeEvent(input$heatMapfileLoaded, {
        HeatMapMatrix1 <<- read.csv(input$heatMapfileLoaded$datapath,header=TRUE,row.names = 1)
        output$heatmap1 <- renderD3heatmap({
          if ( is.null(HeatMapMatrix1))
            return(NULL)
          d3heatmap(
            HeatMapMatrix1,
            colors = input$palette1,
            xaxis_height = 200,
            yaxis_width = 120,
            digits = 4L,
            dendrogram = input$dendrogram1
          )
        })
      })
      
      # Read each file for creating heatmap 
       
      observeEvent(input$fileLoaded, {
        shinyjs::enable("heatmapmatrix")
        shinyjs::disable("saveFile")
        fileName <- input$fileLoaded$name
        fileNameList <<-c(fileNameList, input$fileLoaded$name)
        filePathList<<-c(filePathList, input$fileLoaded$datapath)
        if  (is.null(fileName)){
          fileName <- character(0)
        }  else {
          # Set the files and select items
          if ( grepl("autism_prediction.v2.rscores.exp100",fileName)){
            fileAutismData <<- input$fileLoaded$datapath
          } else {
            if ( grepl("brainspan.spatio-temporal.rma.st-rowmodz-stoufferzge2",fileName)){
              fileBrainData <<- input$fileLoaded$datapath
            } else {
              if ( grepl("reg-per.sub",fileName)){
                fileRegData <<- input$fileLoaded$datapath
              }}}
          if(!is.null(fileAutismData) && !is.null(fileBrainData) && !is.null(fileRegData))
            shinyjs::enable("heatmapmatrix")
          updateSelectInput(session, "inSelect",
                            label = paste("Select Input File"),
                            choices = c(fileName, fileNameList),
                            selected = tail(fileName, 1)
          )
        } 
      })
    } # server
  ))} #runapp

