library(shiny)
library(shinyjs)

source("./utils.R", local = T)

options(shiny.maxRequestSize = 1024^3)  # Maximum upload size 1 GB

server <- function(input, output, session) {
  # >>> configure for deploy >>>
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # ----------------------------------------------------------------
  # -------------------- 根据输入格式更新选项 ----------------------
  # ----------------------------------------------------------------
  # >>> 更新DAM或QH >>>
  output$dam.qh <- renderUI({
    input.format <- input$input.format
    if (input.format == "Proteome Discoverer") {
      choices <- c("DDA")
      selectInput("dam.qh", "Data aquisition method:", choices = choices)
    } else if (input.format == "DIA-NN") {
      choices <- c("DDA", "DIA")
      selectInput("dam.qh", "Data aquisition method:", choices = choices)
    } else if (input.format == "MaxQuant") {
      choices <- c("LFQ in MaxQuant output", "LFQ calculated by iq")
      selectInput("dam.qh", "Quantity by:", choices = choices)
    } else if (input.format == "Spectronuat") {
      choices <- c("FG.MS1Quantity", "FG.MS2Quantity")
      selectInput("dam.qh", "Quantity by:", choices = choices)
    }
  })
  # >>> 更新定量方法输入框 >>>
  output$quantity.method <- renderUI({
    input.format <- input$input.format
    dam.qh <- input$dam.qh
    try({
      if (input.format == "MaxQuant" & dam.qh == "LFQ in MaxQuant output") {
        
      } else {
        selectInput("quantity.method", "Quantity method:", 
                    choices = c("MaxLFQ", "TopN", "Mean", "Median Polish"), 
                    selectize = T)
      }
    }, silent = T)
  })
  # >>> 更新输入文件框 >>>
  output$file.input <- renderUI({
    input.format <- input$input.format
    dam.qh <- input$dam.qh
    try({
      if (input.format == "MaxQuant") {
        card(      
          fileInput("data.protein.group", "Protein group file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv")), 
          fileInput("data.peptide.group", "Peptide group file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv")), 
          fileInput("data.msms", "MS/MS file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv")), 
          fileInput("data.evidence", "Evidence file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv")), 
          fileInput("data.modification", "Oxidation modification file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv"))
        )
      } else if (input.format == "Proteome Discoverer") {
        card(
          fileInput("data.protein", "Protein file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv")), 
          fileInput("data.peptide.group", "Peptide group file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv")), 
          fileInput("data.msms", "MS/MS file",
                    multiple = FALSE,
                    accept = c(".txt", ".tsv"))
        )
      } else {
        fileInput("data.report", "Report (long format)",
                  multiple = FALSE,
                  accept = c(".txt", ".tsv"))
      }
    }, silent = T)
  })
  
  # ----------------------------------------------------------------
  # ------------------------- 执行程序 -----------------------------
  # ----------------------------------------------------------------
  # >>> read data >>>
  
  source("loadData_server.R", local = T)
  report <- reactiveValues(status = FALSE, warning.test = FALSE, suggest.text = "")
  output.report <- eventReactive(params(), {
    if (params()$input.format == "Spectronuat") {
      source("spectronuat_server.R", local = T)
    } else if (params()$input.format == "DIA-NN") {
      source("diann_server.R", local = T)
    } else if (params()$input.format == "Proteome Discoverer") {
      source("pd_server.R", local = T)
    } else if (params()$input.format == "MaxQuant"){
      source("maxquant_server.R", local = T)
    }
    return(output.info())
  })
  
  observeEvent(input$run.button, {print(output.report())})
  
  output.file <- eventReactive(output.report(), {
    if (!is.null(output.report()$warning.test)) {
      output$warning.text <- renderText(output.report()$warning.test)
      output$suggest.text <- renderText(output.report()$suggest.test)
      return(NULL)
    } else {
      return(output.report())
    }
  })
  
  tryCatch({
    output$data <- renderTable({
      head(output.file()$protein.table, 10)
    })
  })

  # >>> download >>>
  source("download_server.R", local = T)
  

}