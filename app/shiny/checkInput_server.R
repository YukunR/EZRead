# sourced by '*_server.R'

# Required input:
## params()
## output

# check input file

#' @param data 各个来源的数据
#' @param dam.qh 数据采集方式或定量头, data aquisition method, character, one of DDA and DIA
#' @param format 数据格式, character, one of pd, mq, diann, sp, autodetect
checkFormat <- function(data, dam.qh, format = "autodetect") {
  if (all(colnames(data)[1: 2] == c("PSMs.Workflow.ID", "PSMs.Peptide.ID"))) {
    res <- list(data = data, dam = dam.qh, format = "Proteome Discoverer")
  } else if (all(colnames(data)[1: 2] == c("Sequence", "Length"))) {
    res <- list(data = data, dam = dam.qh, format = "MaxQuant")
  } else if (all(colnames(data)[1: 2] == c("File.Name", "Run"))) {
    res <- list(data = data, dam = dam.qh, format = "DIA-NN")
  } else if (all(colnames(data)[1: 2] == c("R.Data.Points.per.Peak..MS2..EXT.", "R.MS1.Average.Tolerance..ppm."))) {   # TODO: 标准化所有SP结果
    res <- list(data = data, dam = dam.qh, format = "Spectronuat")
  } else {
    res <- list(error = paste0(
      "Unknown file format with colnames start with ", colnames(data)[1: 2], "."
    ))
  }
  return(res)
}

empty.input <- names(params()[sapply(params(), FUN = is.null)])
if (length(empty.input) > 0) {
  output$warning.text <- renderText(paste0("Input(s) ", paste(empty.input, collapse = ", "), "  are empty!"))
  output$suggest.text <- renderText("Please check your input!")
  req("")
}

# TODO: 添加更多检查
if (params()$input.format == "Proteome Discoverer") {
  format.check <- checkFormat(data.msms, params()$dam.qh, format = params()$input.format)
} else if (params()$input.format == "DIA-NN") {
  format.check <- checkFormat(data.report, params()$dam.qh, format = params()$input.format)
} else if (params()$input.format == "MaxQuant") {
  format.check <- checkFormat(data.evidence, params()$dam.qh, format = params()$input.format)
} else if (params()$input.format == "Spectronuat") {
  format.check <- checkFormat(data.report, params()$dam.qh, format = params()$input.format)
}

if (!is.null(format.check$error)) {
  output$warning.text <- renderText(format.check$error)
  output$suggest.text <- renderText("Please check your input!")
  req("")
}

