# sourced by 'server.R'

# Required input:
## input$download.button


observe({
  if(!is.null(output.file())){
    shinyjs::enable("download.button")
  } else {
    shinyjs::disable("download.button")
  }
})



output$download.button <- downloadHandler(
  filename = function() {
    paste0("output-", params()$input.format, " file quantified by ", params()$quantity.method, "-", Sys.time(), ".xlsx")
  },
  content = function(file) {
    if(report$status == FALSE) NULL else {
      wb_save(output.file()$result.table, file)
    }
    # TODO: 先存下来，再复制
    
    # plot.output.path <- paste0("./output-", input$method, "-", input$metabolism.data.file[1], ".pdf")
    # ggsave(plot.output.path, local.output$plot.res, width = 10, height = 5)
    # data.output.path <- paste0("./output-", input$method, "-", input$metabolism.data.file[1], ".csv")
    # data.write <- data.frame(Compounds = row.names(local.output$data.res), local.output$data.res)
    # write.csv(data.write, data.output.path, quote = T, row.names = F)
    # 
    # zip(zipfile = file, files = c(plot.output.path, data.output.path))
  }
)