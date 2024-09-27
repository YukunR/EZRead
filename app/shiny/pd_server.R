# sourced by 'server.R'

# Required input:
## params()

# operate PD input

output.info <- eventReactive(params(), {
  withProgress(message = "Read data... Please wait...", {
    output.info <- list()
    tryCatch({
      print(params()$data.protein$datapath)
      data.protein <- readData(params()$data.protein$datapath)
      data.peptide.group <- readData(params()$data.peptide.group$datapath)
      data.msms <- readData(params()$data.msms$datapath)
    }, error = function(er) {
      report$warning.text <- er[["message"]]  # renderText(er[["message"]])
      report$suggest.text <- "Please check your input!"  # renderText("Please check your input!")
      return(output.info)
    })
    incProgress(0, message = "Check input format... Please wait...")
    source("./checkInput_server.R", local = T)
    
    incProgress(1 / 4, message = "Make protein list... Please wait...")
    tryCatch({
      protein.list <- makeProteinListPD(data.msms, normalization = params()$normalization)
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
    
    # TODO: format convert and plot
    incProgress(1 / 4, message = "Make result table... Please wait...")
    result.list <- makeResultListPD(protein.table, data.protein, data.peptide.group, data.msms)
    result.table <- groupRows(result.list$result.table, result.list$rows.grouping)
    output.info <- list(protein.table = protein.table, result.table = result.table)
    report$status <- T
    incProgress(1 / 4, message = "Done!")
  })
  print(output.info)
  return(output.info)
})