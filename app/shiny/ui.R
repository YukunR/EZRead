library(shiny)
library(bslib)

# Define UI for data upload app ----
ui <- fluidPage(
  navbarPage("EZRead", 
             tabPanel("Main", 
                      "This is main Page"
                      ), 
             tabPanel("Converter", 
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("input.format", "Input format:", 
                                      choices = c("Proteome Discoverer", "DIA-NN", "MaxQuant", "Spectronuat"), 
                                      selectize = T), 
                          tags$hr(),
                          uiOutput("dam.qh"), 
                          tags$hr(),
                          uiOutput("quantity.method"), 
                          tags$hr(), 
                          uiOutput("file.input"), 
                          tags$hr(), 
                          checkboxInput("normalization", "Median Normalzation"), 
                          tags$hr(),
                          layout_columns(
                            actionButton("run.button", "Run"),
                            downloadButton("download.button", "Download")
                            ), 
                          textOutput("warning.text"), 
                          textOutput("suggest.text")
                          ),
                        mainPanel(
                          tags$h5("only show first 10 rows."), 
                          tags$hr(),
                          tableOutput("data"),
                          plotOutput("plot")
                        )
                      )
             
             ), 
             tabPanel("Contact Us!")
             )
)


# ui <- fluidPage(
#   
#   # App title ----
#   titlePanel("EZRead Format Converter V1.0"),
#   
#   # Sidebar layout with input and output definitions ----
#   sidebarLayout(
#     
#     # Sidebar panel for inputs ----
#     sidebarPanel(
#       selectInput("method", "Method:",
#                   choices = c("SERRF", "LOESS", "QCISRF"), 
#                   selectize = T), 
#       
#       # Horizontal line ----
#       tags$hr(),
#       
#       # Input: Select a file ----
#       fileInput("metabolism.data.file", "Metabolism file (.xlsx)",
#                 multiple = FALSE,
#                 accept = c(".xlsx", "xls")),
#       fileInput("metabolism.order.file", "Metabolism Sample Order (.xlsx)",
#                 multiple = FALSE,
#                 accept = c(".xlsx", "xls")), 
#       tags$hr(), 
#       
#       textInput("prefix", "Prefix to delete in Metabolism Sample Order file", value = "", placeholder = "e.g. MQ2103-013-"), 
#       textInput("suffix", "Suffix to delete in Metabolism Sample Order file", value = ""), 
#       tags$hr(),
#       
#       layout_columns(
#         actionButton("run.button", "Run"), 
#         downloadButton("download.button", "Download"), 
#       ), 
#       textOutput("warning.text"),
#       textOutput("suggest.text")
#     ),
#     
#     # Main panel for displaying outputs ----
#     mainPanel(
#       
#       # Output: Data file ----
#       navset_card_tab(
#         nav_panel("Data", 
#                   tags$h5("only show first 10 rows."), 
#                   tags$hr(), 
#                   tableOutput("data")), 
#         nav_panel("Plot", plotOutput("plot"))
#         
#       )
#     )
#   ))


