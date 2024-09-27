# sourced by 'server.R'

# Required input:
## input$run.button


params <- eventReactive(input$run.button, {
  input.format <- input$input.format
  dam.qh <- input$dam.qh
  quantity.method <- input$quantity.method
  normalization <- input$normalization
  
  param.list <- list("input.format" = input.format, "dam.qh" = dam.qh, 
                     "quantity.method" = quantity.method, "normalization" = normalization)
  
  if (input.format == "Proteome Discoverer") {
    param.list <- c(param.list, list("data.protein" = input$data.protein))
    param.list <- c(param.list, list("data.peptide.group" = input$data.peptide.group))
    param.list <- c(param.list, list("data.msms" = input$data.msms))
  } else if (input.format == "DIA-NN") {
    param.list <- c(param.list, list("data.report" = input$data.report))
  } else if (input.format == "MaxQuant") {
    param.list <- c(param.list, list("data.protein.group" = input$data.protein.group))
    param.list <- c(param.list, list("data.peptide.group" = input$data.peptide.group))
    param.list <- c(param.list, list("data.msms" = input$data.msms))
    param.list <- c(param.list, list("data.evidence" = input$data.evidence))
    param.list <- c(param.list, list("data.modification" = input$data.modification))
  } else if (input.format == "Spectronuat") {
    param.list <- c(param.list, list("data.report" = input$data.report))
  }
  
  return(param.list)
})