
library(shiny)
library(ggplot2)
library(reshape2)
library(DT)

ui <- fluidPage(
  titlePanel("Correlations for annotated neutral masses"),
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_mass", "Selecciona una Masa Neutra:", choices = NULL),
      downloadButton("downloadPlots", "Download Plots")
    ),
    mainPanel(
      plotOutput("heatmap"),
      dataTableOutput("correlation_table")
    )
  )
)

server <- function(input, output, session) {
  correlation_matrices <- correlations_list_2
  updateSelectInput(session, "selected_mass", choices = names(correlation_matrices))
  
  output$heatmap <- renderPlot({
    req(input$selected_mass)
    cor_matrix_2 <- correlation_matrices[[input$selected_mass]]
    cor_data <- as.data.frame(as.table(as.matrix(cor_matrix_2)))
    colnames(cor_data) <- c("Var1", "Var2", "value")
    ggplot(cor_data, aes(Var1, Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "magenta", high = "turquoise", mid = "white", midpoint = 0) +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
      theme_minimal() +
      labs(title = paste("Correlations for Neutral Mass:", input$selected_mass), x = NULL, y = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$correlation_table <- renderDataTable({
    req(input$selected_mass)
    cor_matrix <- correlation_matrices[[input$selected_mass]]
    cor_df <- as.data.frame(cor_matrix)
    cor_df
  }, options = list(pageLength = 10))
  
  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste("plots", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)  # Open the PDF device
      for (mass in names(correlation_matrices)) {
        cor_matrix <- correlation_matrices[[mass]]
        cor_data <- as.data.frame(as.table(as.matrix(cor_matrix_2)))
        colnames(cor_data) <- c("Var1", "Var2", "value")
        
        p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
          geom_tile() +
          scale_fill_gradient2(low = "magenta", high = "turquoise", mid = "white", midpoint = 0) +
          geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
          theme_minimal() +
          labs(title = paste("Correlations for Neutral Mass:", mass), x = NULL, y = NULL) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)  # Print each plot to the PDF
      }
      dev.off()  # Close the PDF device
    }
  )
}

shinyApp(ui = ui, server = server)
