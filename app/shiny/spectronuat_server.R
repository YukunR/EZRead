# sourced by 'server.R'

# Required input:
## params()

# operate SP input

source("./spectronuat.R", local = T)

output.info <- observeEvent(params(), {
  withProgress(message = "Read data... Please wait...", {
    output.info <- list()
    tryCatch({
      if (!is.null(params()$data.report)) {
        data.report <- readData(params()$data.report$datapath)
      }
    }, error = function(er) {
      report$warning.text <- er[["message"]]  # renderText(er[["message"]])
      report$suggest.text <- "Please check your input!"  # renderText("Please check your input!")
      return(output.info)
    })
    
    incProgress(0, message = "Check input format... Please wait...")
    source("checkInput_server.R", local = T)
    
    incProgress(1 / 4, message = "Make protein list... Please wait...")
    tryCatch({
      protein.list <- makeProteinList(data.report, params()$input.format, normalization = params()$normalization, intensity.col = params()$dam.qh)
    }, error = function(er) {
      report$warning.text <- er[["message"]]
      report$suggest.text <- "Please check your input!"
      return(output.info)
    })
    
    incProgress(1 / 4, message = "Make protein table... Please wait...")
    tryCatch({
      protein.table <- createProteinTable(protein.list, method = params()$quantity.method)
    }, error = function(er) {
      report$warning.text <- er[["message"]]
      report$suggest.text <- "Please check your input!"
      return(output.info) 
    })
    incProgress(1 / 4, message = "Make result table... Please wait...")
    result.list <- makeResultListSP(protein.table, data.report)
    groupRows(result.list$result.table, result.list$rows.grouping)
    output.info <- list(protein.table = protein.table, result.table = result.table)
    report$status <- T
    incProgress(1 / 4, message = "Done!")
  })
  return(output.info)
})