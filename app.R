library(shiny)
library(ggplot2)
library(tidyverse)
library(stats)
library(cowplot)
########################
#### User Interface ####
########################

ui <- pageWithSidebar(
  
  # App title ----
  titlePanel(title = div(img(src = "logo.jpg", height = 100, width = 140), "Phage-ELF: Phage Estimator of Lytic Function")),
  
  # Sidebar panel for inputs ----
  sidebarPanel(
    
    fileInput("total_df", "Upload DLS data for model training and testing",
              accept = c(".csv")
    ),
    
    fileInput("size_bins", "Upload size bins for DLS data (Optional)",
              accept = c(".csv")
    ),
    
    fileInput("comp_pairs", "Upload training control-treatment pairs",
              accept = c(".csv")
    ),
    
    fileInput("lytic", "Upload lytic activity data for training",
              accept = c(".csv")
    ),
    
    fileInput("test_pairs", "Upload testing control-treatment pairs (Optional)",
              accept = c(".csv")
    ),
    
    tags$hr(),
    
    selectInput("metric", "Metric used",
                c("Intensity" = "intens",
                  "Volume" = "volume",
                  "Number" = "number")),
    tags$hr(),
    
    strong("Download training data table"),
    br(),
    downloadButton("downloadTrain", "Download CSV file"),
    br(),
    br(),

    strong("Download prediction table"),
    br(),
    downloadButton("downloadTest", "Download CSV file"),
    br()

  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    tabsetPanel(
      tabPanel("Instructions",
               br(),
               p("This webpage accompanies the paper 'Bacteriophage Enumeration and Bioactivity Prediction Using Dynamic Light Scattering'. 
    It provides a user-friendly platform for scientists to train a model with paired DLS and lytic activity data and use this model to
    predict lytic activity from new DLS data."),
               p("The user must upload CSV files with DLS data, lytic data and treatment-control pairs to train the model. CSV files can 
                 be obtained from the 'Save as' option in Excel sheets."),
               br(),
               p("The DLS file must contain a 'Sample.Name' column, followed
                 by columns with Intensity, Volume and Number data. Column names must contain one of these metric names in them for 
                 the webpage to process the data properly (see example below). Additional columns in the DLS output do not need to be 
                 removed."),
               tableOutput('DLS_example'),
               p("An optional file with a single 'size' column can be uploaded to provide the size bins corresponding to each column. This does not affect the model,
               and is only used to scale the x-axis of the DLS figures. If provided, the number of size bins must match the number of metric columns in the DLS data."),
               tableOutput('size_example'),
               p("The model requires a CSV file with a 'control' and 'treatment' columns outlining which samples
               to compare to each other. Sample names must match the names in the DLS data. Include a row where the control is compared to itself for model robustness and as a sanity check."),
               tableOutput('training_example'),
               p("Finally, model training requires lytic data for each control and treatment in the following format:"),
               tableOutput('lytic_example'),
               p("The 'Training data', 'DLS graph' and 'Model' should be populated after lytic data is provided. The user can then provide control-treatment pairs
                 to predict lytic activity from DLS data with the trained model. DLS data for the testing pairs should be included in the DLS csv and the testing pairs should be formatted 
                 like the training pairs:"),
               tableOutput('testing_example'),
               p("If the webpage displays an error message, double-check column names and make sure they match the ones described here. All code for this Shiny app is available
                 at https://github.com/jpourtois/DLS_analysis. In addition, we provide a R software package for more flexibility and higher throughput available for download at https://github.com/jpourtois/phageELF. 
                 Questions can be directed to jp22@alumni.princeton.edu")
               ), 
      tabPanel("Training data",
               br(),
               p("This table summarizes the training data and provides the difference in area-under-the-curve (AUC) for the DLS data 
               and the difference in lytic activity between treatment and control (log10(Treatment) - log10(Control)). For example, a difference of -2 
               represents a decrease in lytic activity of two orders of magnitude. This table can be downloaded in the sidebar."),
               tableOutput("training_data")),
      tabPanel("DLS graph", 
               br(),
               p('These plots compare DLS curves for the treatment and control defined in the control-treatment pairs.'),
               plotOutput("plot")),
      tabPanel("Model",
               br(),
               p('The difference in area-under-the-curve is measured from the DLS curves and plotted against the difference in lytic activity between the treatments
                 and control. A linear model with a fixed intercept at 0,0 is then run. If test DLS data is provided, this linear model can then be used to make predictions
                 about the loss of lytic activity.'),
               plotOutput('linear')),
      tabPanel("Predictions",
               br(),
               p('This table provides the predicted loss of lytic activity if test data is provided. It can be downloaded in the sidebar.'),
               tableOutput("AUC_test_table"))
    )
    
  )
)

################
#### Server ####
################

server <- function(input, output) {
  
  # Read and process all files ----
  
  # Read DLS data and restrict to relevant metric
  data <- reactive({
    
    inFile <- input$total_df
    if (is.null(inFile)) { return(NULL) }    
    dataFile <- read.csv(inFile$datapath)
    
    metric <- input$metric
    
    size_info <- dataFile[,grep(metric,colnames(dataFile), ignore.case = TRUE)]
    
    df_dls <- cbind(dataFile[,c('Sample.Name')],size_info)
    
    colnames(df_dls) <- c('sample',colnames(df_dls)[-1])
    
    df_dls <- df_dls[,1:(ncol(df_dls))]
    
    return(df_dls)
    
  })
  
  # Read size file
  size_df <- reactive({
    
    size_file <- input$size_bins
    
    if (is.null(size_file)) { return(NULL) }    
    size <- read.csv(size_file$datapath)
    
  })
  
  # Read control-treatment pairs for training model
  comp_training <- reactive({
    
    comp_pairs_file <- input$comp_pairs
    
    if (is.null(comp_pairs_file)) { return(NULL) }    
    comp_pairs <- read.csv(comp_pairs_file$datapath)
    
  })
  
  # Read titer data for training model
  titer_data <- reactive({
    
    inFile <- input$lytic
    if (is.null(inFile)) { return(NULL) }    
    dataFile <- read.csv(inFile$datapath)
    return(dataFile)
    
  })
  
  # Read DLS data for model inference
  test_data <- reactive({
    
    inFile <- input$test_data
    if (is.null(inFile)) { return(NULL) }    
    dataFile <- read.csv(inFile$datapath)
    return(dataFile)
    
  })
  
  # Read control-treament pairs for model inference
  test_comp <- reactive({
    
    inFile <- input$test_pairs
    if (is.null(inFile)) { return(NULL) }    
    dataFile <- read.csv(inFile$datapath)
    return(dataFile)
    
  })
  
  # Processing ----
  
  ## Function to calculate AUC difference
  
  compute_AUC <- function(dls_data, comp_pairs, metric){
    
    size <- dls_data[,grep(metric,colnames(dls_data), ignore.case = TRUE)]
    
    diff <- tibble(treatment = comp_pairs$treatment, AUC = NaN)
    
    for (n in 1:nrow(comp_pairs)) {
      
      control <- size[dls_data$sample == comp_pairs$control[n],]
      treatment <- size[dls_data$sample == comp_pairs$treatment[n],]
      diff$AUC[n] <- sum(abs(control-treatment))
    }
    
    return(diff)
    
  }
  
  ## Do analysis on the training data
  
  # Calculate AUC difference 
  AUC <- reactive({
    
    dls_data <- data()
    comp_pairs <- comp_training()
    metric <- input$metric
    
    return(compute_AUC(dls_data,comp_pairs,metric))
    
  })
  
  # Calculate titer difference
  titer_diff <- reactive({
    
    auc_diff <- AUC()
    titer <- titer_data()
    comp_pairs <- comp_training()
    
    titer_dls <- merge(auc_diff,titer, by = 'treatment')
    titer_dls$lytic_diff <- NaN
    
    for (n in 1:nrow(titer_dls)) {
      
      control <- titer$lytic[comp_pairs$control[comp_pairs$treatment == titer_dls$treatment[n]] == titer$treatment]
      treatment <- titer$lytic[titer$treatment == titer_dls$treatment[n]]
      titer_dls$lytic_diff[n] <- log10(treatment+1) - log10(control+1)
    }
    
    return(titer_dls)
    
  })
  
  # Format AUC and titer loss data to display in table
  training_data_table <- reactive({
    
    df <- titer_diff()
    
    df$lytic <- format(df$lytic, scientific = TRUE) 
    
    colnames(df) <- c('Treatment','AUC difference','Titer','Titer loss')
    
    return(df)
    
  })
  
  # Calculate AUC for test data
  AUC_test <- reactive({
    
    dls_data <- data()
    comp_pairs <- test_comp()
    metric <- input$metric
    
    if (is.null(dls_data)|is.null(comp_pairs)){
      output <- NULL
    } else {
      output <- compute_AUC(dls_data,comp_pairs,metric)
    }
    
    return(output)
    
  })
  
  # Predict titer loss from AUC for test data and make table
  test_data_table <- reactive({
    
    titer_dls <- titer_diff()
    
    AUC_test_df <- AUC_test()
    
    linear.model <- lm(formula = lytic_diff ~ AUC - 1, data = titer_dls)
    
    AUC_test_df$pred <- predict(linear.model, AUC_test_df, type='response', se = TRUE)$fit
    AUC_test_df$SE <- predict(linear.model, AUC_test_df, type='response', se = TRUE)$se.fit
    
    colnames(AUC_test_df) <- c('Sample','AUC Difference', 'Titer loss pred. (Log10)', 'SE (Log10)')
    
    return(AUC_test_df)
    
  })
  
  
  # OUTPUTS ----
  
  # Example of DLS data
  output$DLS_example <- renderTable({
    example_DLS <- tibble('Sample.Name' = c('LPS_control','LPS_t1','LPS_t2', 'LPS_t3'))
    example_DLS$'Intensities.1.Percent' <- c(0,0,0,0)
    example_DLS$'Intensities.2.Percent' <- c(0.2,0,0.5, 0.4)
    example_DLS$'Intensities.3.Percent' <- c(0.7,0.5,0.2,0.2)
    example_DLS$'Intensities.4.Percent' <- c(0.1,0.5,0.3,0.4)
    example_DLS$'Intensities.5.Percent' <- c(0,0,0,0)
    return(example_DLS)
  })
  
  # Example of size data
  output$size_example <- renderTable({
    example_size <- tibble('size' = 10^(seq(-2,2)))
    return(example_size)
  })
  
  # Example of control-treatment pair training data
  output$training_example <- renderTable({
    
    example_training <- tibble(control = c('LPS_control','LPS_control','LPS_control'), treatment = c('LPS_control','LPS_t1','LPS_t2'))
    
  })
  
  # Example of lytic data 
  output$lytic_example <- renderTable({
    
    example_lytic <- tibble(treatment = c('LPS_control','LPS_t1','LPS_t2'), lytic = c(10^9,4.2*10^6,8.0*10^7))
    example_lytic$lytic <- format(example_lytic$lytic, scientific = TRUE)
    
    return(example_lytic)
  })
  
  # Example of control-treatment test data
  output$testing_example <- renderTable({
    
    example_testing <- tibble(control = 'LPS_control', treatment = 'LPS_t3')
    
    return(example_testing)
  })
  
  # Render table with AUC and titer loss (training) data
  output$training_data <- renderTable({
    req(data())
    req(comp_training())
    req(titer_data())
    
    training_data_table()
  }) 
  
  # Render table with titer loss predictions 
  output$AUC_test_table <- renderTable({
    req(data())
    req(comp_training())
    req(titer_data())
    req(test_comp())
    
    test_data_table()
  }) 
  
  # Render plot of DLS data
  output$plot<-renderPlot({
    
    req(data())
    req(comp_training())
    req(AUC())
    to_plot <- data()
    
    metric <- input$metric
    
    name.metrics <- c('Intensity','Volume','Number')
    names(name.metrics) <- c('intens','volume','number')
    
    y.axis.name <- name.metrics[metric]
    
    to_plot_long <- gather(to_plot, key = "Size", value = "value", colnames(to_plot[2]):colnames(to_plot)[length(colnames(to_plot))] )
    
    ordered <- to_plot_long[order(to_plot_long$sample),]
    
    # Use size data if provided, otherwise just use integers
    if (!is.null(size_df())) {
      
      ordered$xscale <- size_df()$size
    } else {ordered$xscale <- rep(10^((1:(ncol(to_plot) - 1))/ncol(to_plot)), nrow(to_plot))}
    
    # Read treatment-control pairs
    comp_pairs_file <- input$comp_pairs
    
    if (is.null(comp_pairs_file)) { return(NULL) }    
    comp_pairs <- read.csv(comp_pairs_file$datapath)
    
    ordered <- ordered[ordered$sample %in% comp_pairs$treatment,]
    
    controls_name <- unique(comp_pairs$control)
    
    # Make plot for each control and its treatments DLS data
    plot_data_column = function (data, comp_pairs, control) {
      
      ordered <- data
      pairs_to_plot <- comp_pairs[comp_pairs$control == control,]
      to_plot <- ordered[ordered$sample %in% pairs_to_plot$treatment,]
      
      ggplot(to_plot,aes(x=xscale,y=value, color = sample)) +
        geom_line() +
        labs(x = 'Size', y = y.axis.name, color = 'Treatment') +
        scale_x_continuous(trans='log10') +
        theme_bw()
    }
    
    # Make a plot for each control
    myplots <- lapply(controls_name, plot_data_column, data = ordered, comp_pairs = comp_pairs)
    
    plot_grid(plotlist=myplots, ncol = 1)},height = 550,width = 550)
  
  # Render plot with training data, model and test data and predictions if provided
  output$linear <- renderPlot({
    
    req(data())
    req(comp_training())
    req(titer_data())
    
    titer_dls <- titer_diff()
    
    # Train model
    lm.model <- lm(formula = lytic_diff ~ AUC - 1, data = titer_dls)
    
    newdata <- data.frame(AUC = seq(0, 200,len=500))
    
    #use fitted model to predict values of vs
    newdata$pred <- predict(lm.model, newdata, interval = 'confidence')[,'fit']
    newdata$upper <- predict(lm.model, newdata, interval = 'confidence')[,'upr']
    newdata$lower <- predict(lm.model, newdata, interval = 'confidence')[,'lwr']
    
    # Add test data on graph if test data is present
    
    test_tibble = tibble(sample = factor(), AUC = numeric(), pred = numeric())
    
    if (!is.null(AUC_test())){
      
      AUC_test_data <- AUC_test()
      AUC_test_data$pred <- predict(lm.model, AUC_test_data)
      test_tibble <- rbind(test_tibble, AUC_test_data)
      
      ggplot() +
        geom_point(data = titer_dls, aes(x = AUC, y = lytic_diff, color = "Training data")) +
        geom_line(data = newdata, aes(x = AUC, y = pred)) +
        geom_point(data = test_tibble, aes(x = AUC, y = pred, color = "Test data"), shape = 19, size = 2.5)+
        geom_ribbon(data = newdata, aes(x = AUC,ymin = lower, ymax = upper), alpha=0.3) +
        labs(x = 'AUC Difference', y = 'Difference in lytic activity (log10)') +
        theme_bw()
    } else{
      
      ggplot() +
        geom_point(data = titer_dls, aes(x = AUC, y = lytic_diff)) +
        geom_line(data = newdata, aes(x = AUC, y = pred)) +
        geom_ribbon(data = newdata, aes(x = AUC,ymin = lower, ymax = upper), alpha=0.3) +
        labs(x = 'AUC Difference', y = 'Difference in lytic activity') +
        theme_bw()
    }
    
    
  },height = 300,width = 500)
  
  # Allow to download data table
  output$downloadTrain <- downloadHandler(
    
    filename = 'AUC_titer_training.csv',
    content = function(file) {
      write.csv(training_data_table(), file, row.names = FALSE)
    }
  )
  
  # Allow to download predictions
  output$downloadTest <- downloadHandler(
    
    filename = 'AUC_titer_test.csv',
    content = function(file) {
      write.csv(test_data_table(), file, row.names = FALSE)
    }
  )
    
  
}

shinyApp(ui, server)