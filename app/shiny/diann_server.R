# sourced by 'server.R'

# Required input:
## params()

# operate DIANN input

source("./diann.R", local = T)

output.info <- eventReactive(params(), {
  withProgress(message = "Read data... Please wait...", {
    output.info <- list()
    tryCatch({
      if (!is.null(params()$data.report)) {
        print(2)
        data.report <- readData(params()$data.report$datapath)
      }
    }, error = function(er) {
      report$warning.text <- er[["message"]]  # renderText(er[["message"]])
      report$suggest.text <- "Please check your input!"  # renderText("Please check your input!")
      return(output.info)
    })
    
    incProgress(0, message = "Check input format... Please wait...")
    source("checkInput_server.R", local = T)
    
    tryCatch({
      if (params()$quantity.method == "MaxLFQ") {
        incProgress(1 / 2, message = "Quantification... Please wait...")
        protein.table <- diannLFQ(data.report, params()$dam.qh)
      } else {
        incProgress(1 / 4, message = "Make protein list... Please wait...")
        protein.list <- makeProteinList(data.report, params()$input.format, normalization = params()$normalization)
        incProgress(1 / 4, message = "Make protein table... Please wait...")
        protein.table <- createProteinTable(protein.list, method = params()$quantity.method)
      }
    }, error = function(er) {
      report$warning.text <- er[["message"]]  # renderText(er[["message"]])
      report$suggest.text <- "Please check your input!"  # renderText("Please check your input!")
      return(output.info)
    })

    print(params()$input.format)
    incProgress(1 / 4, message = "Make result table... Please wait...")
    result.list <- makeResultListDIANN(protein.table, data.report)
    result.table <- groupRows(result.list$result.table, result.list$rows.grouping)
    output.info <- list(protein.table = protein.table, result.table = result.table)
    report$status <- T
    incProgress(1 / 4, message = "Done!")
  })
  return(output.info)
})