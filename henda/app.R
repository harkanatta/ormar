# Load necessary libraries
library(shiny)
library(DT)
library(readr)

# Load the data (update the path if needed)
KolgrTaxa <- read_csv("/home/rstudio/data/raw/KolgrTaxa.csv", na = "empty")

# Handle invalid UTF-8 in KolgrTaxa columns
KolgrTaxa[] <- lapply(KolgrTaxa, function(x) iconv(x, from = "UTF-8", to = "UTF-8", sub = NA))
KolgrTaxa <- KolgrTaxa[complete.cases(KolgrTaxa), ]

# Define UI
ui <- fluidPage(
  titlePanel("KolgrTaxa Data Filter"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        inputId = "remove_list",
        label = "Enter terms to remove (autocomplete available):",
        choices = NULL,  # No predefined choices; user can type freely
        multiple = TRUE,
        options = list(create = TRUE)
      ),
      checkboxGroupInput(
        inputId = "search_columns",
        label = "Select columns to search in:",
        choices = names(KolgrTaxa),  # Allow user to select any column for search
        selected = "gamalt"  # Default column to search
      ),
      checkboxGroupInput(
        inputId = "selected_columns",
        label = "Select columns to display:",
        choices = names(KolgrTaxa),  # Choices are all column names
        selected = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "gamalt")  # Default display columns
      )
    ),
    mainPanel(
      dataTableOutput("filtered_data")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  output$filtered_data <- renderDataTable({
    # Start with the full data frame
    filtered_data <- KolgrTaxa
    
    # Check if remove_list has input and user selected columns to search in
    if (!is.null(input$remove_list) && length(input$remove_list) > 0 && !is.null(input$search_columns)) {
      # Create a pattern based on the remove_list
      remove_pattern <- paste(input$remove_list, collapse = "|")
      
      # Apply the filter across selected columns
      filtered_data <- filtered_data[!apply(filtered_data[input$search_columns], 1, function(row) {
        any(grepl(remove_pattern, row, ignore.case = TRUE))
      }), ]
    }
    
    # Select columns chosen by the user
    if (!is.null(input$selected_columns)) {
      filtered_data <- filtered_data[, input$selected_columns, drop = FALSE]
    }
    
    # Render the data table
    datatable(filtered_data, options = list(pageLength = 10))
  })
}

# Run the Shiny App
shinyApp(ui = ui, server = server)
