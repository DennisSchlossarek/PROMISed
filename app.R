# Welcome to PROMISed: A PROtein Metabolite Interaction using Size-Separation Data analysis tool

# List of Required packages:
require(shiny)        # For the App to run
require(shinyBS)      # Shiny meeting bootstrap? Allows hovering Tooltips 
require(shinycssloaders)  #
require(shinybusy)    #
require(DT)           # To render interactive Data-Tables
require(stringr)      # Data handling and organization
require(stats)        # General Statistics
require(plyr)         # Data handling and organization
require(MESS)         # Used for MESS::auc (area under the curve) 
require(dplyr)        # 
require(RColorBrewer) # Nice Color pallettes
require(ggsci)        # pal_startrek color-pallette
require(ggplot2)      # For proper Plotting
require(multcompView) # For plotting (here, significance-indices in boxplot)
require(gridExtra)    # Organisation of ggplots
require(grid)         # Organisation of ggplots
require(igraph)       # For Network creation and analysis
require(visNetwork)   # AWESOME interactive network
require(pheatmap)     # AWESOME pretty heatmaps
require(tidyr)
require(pastecs)      # for local maxima
require(zip)
require(eulerr)       # for venn diagrams
require(VennDiagram)

# Load functions from outside R-scripts
source("./Functions_PROMISEed_v.1.0.0.R")
source("./publication_list_HTML.R")

# Load Example Data Files
metabolite_demo_data  <- read.delim("./dummy_metabolite_data.txt", row.names = 1)
protein_demo_data     <- read.delim("./dummy_protein_data.txt", row.names = 1)

# Create Graphical User Interface
ui <- fluidPage( # Start of USER-INTERFACE
  
  br(""),
  
  tags$head(
    tags$style(HTML('#select_previous{background-color:lightgrey}')),
    tags$style(HTML('#select_next{background-color:lightgrey}'))
    
    ),
  
  # Load a rotating SEC-column gif as loading icon
  add_busy_gif(
    timeout = 1000,
    src = "PROMIS rotating.gif",
    height = 80, width = 80,
    position = "full-page",
    overlay_css = "z-index: 1000"
  ),
  
  fluidRow( # UI Header, containing Dataset and row-selectors
    
    column(3, 
           
           # Title Image: PROMIS Processing Logo,           
           tags$a(imageOutput(outputId = "my_image_processing",
                       width = "340%", 
                       height = "157px")),
           helpText("PROMISed review Version: 1.0.2")
           
           ),
    
    column(1,
           
           br(""),
           
           radioButtons(inputId = "dataset_selector", 
                        label = "Dataset displayed",
                        choices = c("Dataset A", "Dataset B"),
                        selected = "Dataset A")
           ),
    
    column(3, 
         
         uiOutput(outputId = "selector_a_ui"),
         uiOutput(outputId = "selector_b_ui")
         
         ),
    
    column(1,
          br(),
          br(),
          uiOutput(outputId = "select_previous_ui"),
          
          uiOutput(outputId = "select_next_ui")
          ),
     
    column(2, offset = 1,
       
         tags$a(imageOutput(outputId = "logo_mpimp", 
                            width = "212%", 
                            height = "100px"),
                target= "_blank")
       
),
column(1,
       
       tags$a(imageOutput(outputId = "logo_minerva", 
                          width = "212%", 
                          height = "100px"),
              target= "_blank")
       
)

           
  ), # End of Header
  
  tabsetPanel(id = "FirstOrder", # Create Outer-most Tab-Structure
      
              tabPanel("Data Upload",           # Tab to Upload data and confirm columns
                       
                       sidebarLayout( # Create Content in a Sidebarlayout: Inputs in SidebarPanel, Output (Plots) in MainPanel
                         sidebarPanel(width = 3,
                                      
                                      # Data-Upload and Column-Selection for Dataset A
                                      wellPanel( 
                                        fileInput(inputId = "file_a",
                                                  multiple = FALSE, 
                                                  label = "Upload Dataset A",
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv")),
                                        
                                        splitLayout(
                                          uiOutput(outputId = "lower_bound_colnames_ui_a"),
                                          
                                          wellPanel(textOutput(outputId = "lower_bound_colnames_view_a"))
                                        ), 
                                        
                                        splitLayout(        
                                          uiOutput(outputId = "upper_bound_colnames_ui_a"),
                                          
                                          wellPanel(textOutput(outputId = "upper_bound_colnames_view_a"))
                                        )
                                      ),
                                      
                                      # Data-Upload and Column-Selection for Dataset B
                                      wellPanel( 
                                        fileInput(inputId = "file_b",
                                                  multiple = FALSE, 
                                                  label = "Upload Dataset B",
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv")),
                                        
                                        splitLayout(
                                          uiOutput(outputId = "lower_bound_colnames_ui_b"),
                                          
                                          wellPanel(textOutput(outputId = "lower_bound_colnames_view_b"))
                                        ), 
                                        
                                        splitLayout(        
                                          uiOutput(outputId = "upper_bound_colnames_ui_b"),
                                          
                                          wellPanel(textOutput(outputId = "upper_bound_colnames_view_b"))
                                        )
                                      ),
                                      
                                      # Button to confirm data and move to next tab
                                      actionButton(inputId = "upload_button", 
                                                   label = strong("Confirm Data Upload"), 
                                                   width = "100%"),
                                      
                                      br(""),
                                      
                                      wellPanel(
                                        uiOutput("download_demo_data_ui")
                                      )
                                      
                         ), #End of SidebarPanel "Data Upload"
                         
                         # Start of mainPanel (This is where the Plots live)
                         mainPanel(
                           
                           fluidRow(column(8, offset = 1,
                                           br("")
                           )
                           )
                           
                         ) # End of mainPanel  "Data Upload"
                         
                       ) # End of SidebarLayout  "Data Upload"
 
              ), # End of Data-Upload Tab
  
              tabPanel("1. Pre-Processing",
                       
                       sidebarLayout(
                         sidebarPanel(width = 3, #style = "background-color: #28A9A0",
                                      
                                      h4("Pre-Processing"),
                                      
                                      wellPanel(
                                        checkboxInput(inputId = "box_filterpeak", 
                                                      label = "Remove Single Peaks",
                                                      value = TRUE),
                                        
                                        bsPopover(id = "box_filterpeak",
                                                  title = "",
                                                  content = "Remove Peaks Spanning Only a Single Fraction", 
                                                  placement = "right", 
                                                  trigger = "hover"),
                                        
                                        hr(),
                                        
                                        checkboxInput(inputId = "box_max", 
                                                      label = "Normalize",
                                                      value = TRUE),
                                        
                                        bsPopover(id = "box_max", 
                                                  title = "",
                                                  content = "Normalize Fractionation Profiles to Maximum Intensity", 
                                                  placement = "right", 
                                                  trigger = "hover"),
                                        
                                        numericInput(inputId = "num_lower", 
                                                     label = "Minimum Relative Intensity",
                                                     value = 0.1,
                                                     min = 0,
                                                     max = 1, 
                                                     step = 0.05,
                                                     width = "200px"),
                                        
                                        bsPopover(id = "num_lower", 
                                                  title = "",
                                                  content = "Apply a Minimum Relative Intensity Threshold to Remove Low-Intensity-Data-Noise", 
                                                  placement = "right", 
                                                  trigger = "hover"),
                                        
                                        hr(),
                                        
                                        checkboxInput(inputId = "box_loess", 
                                                      label = "Profile Smoothing",
                                                      value = TRUE),
                                        
                                        bsPopover(id = "box_loess",
                                                  title = "",
                                                  content = "Smooth Profiles by Fiting a Polynomial Surface Determined by Span-Value Using Local Fitting", 
                                                  placement = "right", 
                                                  trigger = "hover"),
                                        
                                        numericInput(inputId = "num_span", 
                                                     label = "Span-Value",
                                                     value = 0.15, 
                                                     min = 0.1, 
                                                     max = 1, 
                                                     step = 0.05,
                                                     width = "200px"),
                                        
                                        bsPopover(id = "num_span", 
                                                  title = "",
                                                  content = "Relative Size of Local Neghborhood to Include in Smoothing", 
                                                  placement = "right", 
                                                  trigger = "hover"),
                                        
                                        br(),
                                        
                                        actionButton(inputId = "button_treat", 
                                                     label = strong("Start Analysis"),
                                                     width = "100%")
                                        
                                      ),
                                      
                                      br(),
                                      
                                      uiOutput(outputId = "down_smooth_ui")
                                      
                         ), # End of SidebarPanel 1
                         mainPanel(width = 9,
                                   
                                   tabsetPanel(
                                     tabPanel("Raw Profiles",
                                              br(),
                                              fluidRow(
                                                column(8, offset = 1,
                                                       plotOutput(outputId = "plot", 
                                                                  height = "500px", 
                                                                  width = "100%")
                                                ),
                                                column(2,  offest = 1,
                                                       br(""),
                                                       uiOutput(outputId = "line_plot_check_ui")
                                                )
                                              )
                                     ) 
                                   )
                         ) # End of MainPanel 1
                       )
                       
              ), # End of Tab "Pre-Processing"
              
              tabPanel("2. Replicate Pooling",
                       
                       sidebarLayout(
                         sidebarPanel(width = 3, 
                                      
                                      h4("Replicate Pooling"),
                                      
                                      wellPanel(width = 12,
                                                
                                                radioButtons(inputId = "cor_method_replicate", 
                                                             label = "Correlation Method",
                                                             choiceNames = c("Pearson", "Kendall's-Tau", "Spearman Rank"),
                                                             choiceValues = c("pearson", "kendall", "spearman"),
                                                             selected = "pearson", 
                                                             inline = TRUE
                                                ),
                                                
                                                sliderInput(inputId = "pcc_thresh_replicate", 
                                                            label = "Reproducibility Threshold",
                                                            min = 0, 
                                                            max = 0.99, 
                                                            value = 0.7, 
                                                            step = 0.05),
                                                
                                                bsPopover(id = "pcc_thresh_replicate", 
                                                          title = "", 
                                                          content = "Correlation-Coefficient Threshold to Consider Replicates as Reproducible",
                                                          placement = "top",
                                                          trigger = "hover"),
                                                
                                                checkboxInput(inputId = "box_rep", label = "Keep Single Replicates?",
                                                              value = FALSE),
                                                
                                                bsPopover(id = "box_rep", 
                                                          title = "", 
                                                          content = "Should Entries Only Measured in One Replicate be Kept for Further Analysis?",
                                                          placement = "right",
                                                          trigger = "hover")
                                                ),
                                      
                                      br(),
                                      uiOutput(outputId = "down_rep_comb_ui")
                                      ), 
                         mainPanel(width = 9,
                                   
                                   tabsetPanel(
                                     tabPanel("Pooled Replicates",
                                              
                                              br(),
                                              fluidRow(
                                                column(8, offset = 1,
                                                       plotOutput(outputId = "repcomb_plot", 
                                                                  height = "500px", 
                                                                  width = "100%")
                                                )
                                              )
                                     ),
                                     tabPanel("Global Fractionation Profile",
                                              br(""),
                                              fluidRow(
                                                column(6, offset = 1,
                                                       plotOutput(outputId = "heatmap_plot",
                                                                  height = "600px",
                                                                  width = "100%")
                                                )
                                              )
                                     ) # End of "Global Elution Profile"
                                   ) # End of RepComb/mainPanel/Tabsetpanel
                         ) # End of Repcomb/mainpanel
                       ) # End of Repcomb/Siderbarlayout
              ), # End of RepComb Panel
              
              tabPanel("3. Deconvolution",
                       
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      
                                      h4("Deconvolotuion Settings"),
                                      
                                      wellPanel(width = 12,
                                                
                                               #p(strong("Deconvolotuion Settings")),
                                                
                                                sliderInput(inputId = "var_limit_small", 
                                                            label = "Minimum Relative Intensity",
                                                            min = 0, 
                                                            max = 0.99, 
                                                            value = 0.2, 
                                                            step = 0.05),
                                                
                                                bsPopover(id = "var_limit_small", 
                                                          title = "",
                                                          content = "Apply a Minimum Relative Intensity Threshold to Remove Low-Intensity-Data-Noise", 
                                                          placement = "top", 
                                                          trigger = "hover"),
                                                
                                                sliderInput(inputId = "var_min_peak", 
                                                            label = "Minimum Intensity of Local Maxima", 
                                                            min = 0, 
                                                            max = 0.99, 
                                                            value = 0.2, 
                                                            step = 0.05),
                                                
                                                bsPopover(id = "var_min_peak", 
                                                          title = "",
                                                          content = "Relative Intensity to Consider Local Maxima as Independent Peaks", 
                                                          placement = "top", 
                                                          trigger = "hover"),
                                                
                                                sliderInput(inputId = "var_limit_large", 
                                                            label = "Minimum Incline", 
                                                            min = 0,
                                                            max = 0.99, 
                                                            value = 0.8, 
                                                            step = 0.05),
                                                
                                                bsPopover(id = "var_limit_large", 
                                                          title = "",
                                                          content = "Minimum Increase or Decrease required to Start a New Peak", 
                                                          placement = "top", 
                                                          trigger = "hover")
                                      ),
                                      
                                      checkboxInput(inputId = "nodeconvolution_input", label = "No Deconvolution"),
                                      
                                      br(),
                                      
                                      downloadButton(outputId = "down_decon", 
                                                     label = "Download Deconvoluted Data",
                                                     style = "width:100%;")
                                      
                         ), # End of SidebarPanel 3
                         mainPanel(width = 9,
                                   
                                   tabsetPanel(id = "Deconvo_Tabs",
                                               
                                               tabPanel("Deconvolution")
                                               
                                   ) #End of Deconvolution/Mainpanl/Tabsetpanel
                         ) # End of Deconvolution/MainPanel
                       ) # End of Deconvolution/Siderbarlayout
              ), # End of Deconvolution
              
              tabPanel("4. Data Integration", 
                       
              sidebarLayout(
                sidebarPanel(width = 3,
                             
                             h4("Correlation"),    
                             wellPanel(
                               
                               radioButtons(inputId = "cor_method_integration", 
                                            label = "Correlation Method",
                                            choiceNames = c("Pearson", "Kendall-tau", "Spearman"),
                                            choiceValues = c("pearson", "kendall", "spearman"),
                                            selected = "pearson", 
                                            inline = TRUE
                               ),
                               
                               bsPopover(id = "cor_method_integration",
                                         title = "",
                                         content = "Choose a Method of Correlation for Data Integration", 
                                         placement = "top", 
                                         trigger = "hover"),
                               
                               sliderInput(inputId = "pcc_slider_integration", 
                                           label = "Correlation Threshold", 
                                           value = 0.7, 
                                           min = 0, 
                                           max = 0.99), 
                               
                               bsPopover(id = "pcc_slider_integration",
                                         title = "",
                                         content = "Determine the Threshold of Correlation, <br> Used for Intersection and Network-Mapping", 
                                         placement = "top", 
                                         trigger = "hover"),
                             
                             actionButton(inputId = "button_integration", 
                                          label = "Confirm Correlation Settings", 
                                          width = "100%")
                             ),
                             
                             radioButtons(inputId="euler_input", 
                                           label = "Layout", 
                                           choices  = c("Venn Diagram","Euler Diagram"),
                                           selected = "Venn Diagram",
                                           inline = TRUE),
                             
                             splitLayout(
                              downloadButton(outputId = "down_corr_table_zip", label = "Correlation Tables"),
                             
                              downloadButton(outputId = "down_overlap", label = "Intersection Table")
                             ),
                             
                             bsPopover(id = "down_corr_table_zip",
                                       title = "",
                                       content = "Download Correlation-Tables of All Conditions", 
                                       placement = "right", 
                                       trigger = "hover")
                             
                ), # End of Data Integration/Sidebarlayout/sidebarPanel 
                
                mainPanel(width = 9, # Start of Data Integration/Sidebarlayout/mainPanel
                         
                          tabsetPanel(id = "Integration_Tabs",
                                      
                                      tabPanel("Intersections of Conditions", 
                                               
                                               fluidRow(
                                                 column(5, offset = 1,
                                                        
                                                        br(""), 
                                                        plotOutput(outputId = "venn_set",
                                                                   height = "600px", 
                                                                   width = "600px")
                                                        
                                                       
                                                        
                                                 )
                                                 )
                                               
                                               ) # End of Data Integration/Sidebarlayout/mainPanel/Integration Tabs/Intersection of Conditions
                                      
                                      ) # End of Data Integration/Sidebarlayout/mainPanel/Integration Tabs
                          
                ) # End of Data Integration/Sidebarlayout/mainPanel
              )  # End of Data Integration/Sidebarlayout
              ), # End of Data Integration
              
              tabPanel("5. Network Analysis", value = 5,
                       
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      
                                      div(
                                        h4("Network"),
                                        
                                        wellPanel(
                                          
                                          radioButtons(inputId = "network_filter",
                                                       label = "Filter Network on",
                                                       choiceNames = c("No Filter", "Selection 1", "Selection 2"),
                                                       choiceValues = c("network_filter_0", "network_filter_1",  "network_filter_2"),
                                                       selected = "network_filter_0",
                                                       inline = TRUE),
                                          
                                          hr(),
                                          
                                          splitLayout(
                                            radioButtons(inputId = "node_color_input",
                                                         label = "Node Colors",
                                                          choiceNames = c("Cluster", "k-Coreness"), 
                                                          choiceValues = c(1,2),
                                                         inline = FALSE),
                                            
                                            
                                            radioButtons(inputId = "my_layout", 
                                                         label = "Layout", 
                                                         choiceNames = c("Force-directed","Circles"),
                                                         choiceValues =c(1:2),
                                                         inline = FALSE)
                                            
                                          )
                                        )
                                      )
                                      
                                      ), # End of Network Analysis/Sidebar
                         
                         mainPanel(width = 9,
                                   
                                   tabsetPanel(id = "Network_Tabs",
                                               
                                               tabPanel("Networks")
                                   )
                                   
                                   ) # End of Network Analysis/mainPanel
                       ) # End of Network Analysis/sidebarLayout
              ), # End of Network Analysis
              
              tabPanel("6. Differential Fractionation",
                       
                       sidebarLayout( # Siderbar-Layout Dis-Elution-Score
                         sidebarPanel(width = 3,
                                   
                                               h4("Differntial Fractionation"),
                                               
                                               wellPanel(
                                                 numericInput(inputId = "pvalue", 
                                                              label = "p-Value Threshold",
                                                              value = 0.05, min = 0.001, max = 1, step = 0.005),
                                                
                                                 downloadButton(outputId = "down_table", 
                                                                label = "Download Dis-Elution-Score Results", 
                                                                style = "width:100%;")
                                               )

                                      
                         ), # End of Dis-Elution-Score/Sidebar
                         mainPanel(width = 9,
                                   
                                   tabsetPanel(id = "Manhattan_Tabs",
                                               tabPanel("Combined Plot")
                                   )
                                   
                         ) # End of  Dis-Elution-Score/MainPanel 2
                       )  # End of Dis-Elution-Score/Sidebarlayout
                       
              ), # End of Differential Elution
              
              tabPanel("About",
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      
                                      author_info_HTML()
                                      
                         ), mainPanel(width = 9, offest = 1,
                                      
                                      publication_list_HTML()
                         )
                       )
              ), # End of "About" Tab
              
             tabPanel("Disclaimer",
                       disclaimer_HTML()
                       )
              
  ) # End of TabsetPanel "First Order"
) # End of UI
####

server <- function(input, output, session){ # Start of SERVER-Mainframe
  
  dir.create("promised_tmp")
  
  promised_tmp_wd <- "./promised_tmp"
  code <- paste(paste(sample(x = LETTERS, 2), collapse = ""), sample(999, 1), sep = "")
  
  #setwd(promised_tmp_wd)
  #message_log <- file(paste0(Sys.Date(), "_", code,"_log.txt"), open = "a")
  #sink(file = message_log, type = "message")
  
  observeEvent({input$file_a
                input$file_b}, {   # Observe Uploaded Files A
    
    my_data_upload_a <- data.frame(read.delim(input$file_a$datapath, 
                                            row.names = 1, 
                                            check.names = FALSE))
    
    # Render UI to Select first and last column to be processed:
    output$lower_bound_colnames_ui_a <- renderUI({
      
      div(numericInput(inputId = "lower_bound_colnames_a", 
                       label = "First Column to Analyse", 
                       min = 1,
                       max = ncol(my_data_upload_a), 
                       value = 1))
    })
    
    output$upper_bound_colnames_ui_a <- renderUI({
      
      div(numericInput(inputId = "upper_bound_colnames_a", 
                       label = "Last Column to Analyse", 
                       min = 1, 
                       max = ncol(my_data_upload_a), 
                       value = ncol(my_data_upload_a)))
      
    })
    
    # Render TextOutput for selected colnames for double-checking
  
    output$lower_bound_colnames_view_a <- renderText(colnames(my_data_upload_a)[as.numeric(input$lower_bound_colnames_a)])
    output$upper_bound_colnames_view_a <- renderText(colnames(my_data_upload_a)[as.numeric(input$upper_bound_colnames_a)])
    
    my_data_upload_b <- data.frame(read.delim(input$file_b$datapath, 
                                              row.names = 1, 
                                              check.names = FALSE))
    
    # Render UI to Select first and last column to be processed:
    output$lower_bound_colnames_ui_b <- renderUI({
      
      div(numericInput(inputId = "lower_bound_colnames_b", 
                       label = "First Column to Analyse", 
                       min = 1,
                       max = ncol(my_data_upload_b), 
                       value = 1))
    })
    
    output$upper_bound_colnames_ui_b <- renderUI({
      
      div(numericInput(inputId = "upper_bound_colnames_b", 
                       label = "Last Column to Analyse", 
                       min = 1, 
                       max = ncol(my_data_upload_b), 
                       value = ncol(my_data_upload_b)))
      
    })
    
    # Render TextOutput for selected colnames for double-checking
    
    output$lower_bound_colnames_view_b <- renderText(colnames(my_data_upload_b)[as.numeric(input$lower_bound_colnames_b)])
    output$upper_bound_colnames_view_b <- renderText(colnames(my_data_upload_b)[as.numeric(input$upper_bound_colnames_b)])
   
  observeEvent({input$upload_button}, { # Observe Uploaded Datasets and chosen columns
    
    # Subset both Datasets using the chosen Columns:
    
    my_data_a <- my_data_upload_a[,c(as.numeric(input$lower_bound_colnames_a) : as.numeric(input$upper_bound_colnames_a))]
    my_data_b <- my_data_upload_b[,c(as.numeric(input$lower_bound_colnames_b) : as.numeric(input$upper_bound_colnames_b))]
    
    # This is where it gets messy: 
    # Extract Meta-data (Replicates, Treatments) for both datasets separatly to use for array-dimensions:
    
    metadata_a <- colnames2metadata(my_data_a, input$file_a$name )
    metadata_b <- colnames2metadata(my_data_b, input$file_b$name )
    
    # Render Meta-Data dependent UI Elements:
    
    j <- metadata_a[[2]]
    
    while(j >= 1){
      
      insertTab(inputId = "Deconvo_Tabs",
                
                tabPanel(paste0(metadata_a[[1]][j]), value = j, 
                         
                         br(),
                         
                         fluidRow(
                           column(8, offset = 1,
                                  
                                  plotOutput(outputId = paste0("decon_plot_", j), 
                                             height = "500px", 
                                             width = "100%")
                           )
                         )
                         
                ), target = "Deconvolution", position = "after")
      j <- j - 1  
    }
    
    removeTab(inputId = "Deconvo_Tabs", target = "Deconvolution")
    
    # Insert Tabs for pair-wise dis-elution-score: 
    if(as.numeric( metadata_a[[2]] > 1)){
      
      number_panel <- reactive(ncol(combn(as.numeric( metadata_a[[2]]), 2)))
      name_panel <- namesCombn(metadata_a[[1]])
      
      i <-  number_panel()
      
      while(i >= 1){
        
        insertTab(inputId = "Manhattan_Tabs",
                  
                  tabPanel(paste0(name_panel[i]), value = i, 
                           
                           br(),
                           
                           fluidRow(
                             column(10, offset = 0,
                                    
                                    plotOutput(outputId = paste0("boxplot", i),
                                               height = "700px")
                                    
                             )))
                  , target = "Combined Plot", position = "after")
        
        i <- i - 1  
      }
      removeTab(inputId = "Manhattan_Tabs", target = "Combined Plot")
    }
    
    removeTab(inputId = "FirstOrder", target = "Data Upload")
    
    # Insert Tabs for Data Integration Output:
    # Insert Tabs, dynamicly for number of Treatments:  
    
    k <- metadata_a[[2]]
    
    while(k >= 1){
      
      insertTab(inputId = "Integration_Tabs",
                
                tabPanel(paste0(metadata_a[[1]][k]), value = k, 
                         
                         tabsetPanel(id = "InnerPanels",
                                     
                                     tabPanel("Fractionation Profiles",
                                              br(""),
                                              fluidRow(column(10 ,offset = 1,
                                                              
                                                              plotOutput(outputId = paste0("target_plot_", k),  height = "500px", width = "100%")
                                                              
                                              ))
                                     ), # End of Select-Plot Panel
                                     
                                     tabPanel("Co-Fractionation",
                                              br(""),
                                              fluidRow(column(10 ,offset = 1,
                                                              
                                                              plotOutput(outputId = paste0("hits_plot_", k), height = "500px", width = "100%"),
                                                              
                                                              DTOutput(outputId = paste0("pcc_table_", k))
                                                              
                                              ))
                                     ) # End of Best-Hits Panel
                                     
                         ) # End of Tabsetpanel
                         
                ), # End of Inserted Tabpanel 
                
                
                target = "Intersections of Conditions", position = "after") # End InsertTab
      
      k <- k - 1  
    }
    
    i <- metadata_a[[2]]
    
    while(i >= 1){
      
      insertTab(inputId = "Network_Tabs",
                
                tabPanel(paste0(metadata_a[[1]][i]), value = i, 
                         
                         fluidRow(
                           column(10, offset = 1,
                                  
                                  visNetworkOutput(outputId = paste0("visnetwork_network_", i), 
                                                   height = "600px", 
                                                   width = "900px"),
                                  hr(""),
                                  br(""),
                                  downloadButton(outputId = paste0("network_files_", i),
                                                 label = "Get Network Files"),
                                  
                                  bsPopover(id =  paste0("network_files_", i),
                                            title = "",
                                            content = "Download Edge- and Nodelists as .txt", 
                                            placement = "top", 
                                            trigger = "hover")
                                  
                           )
                         )
                         
                ), # End of Inserted Tabpanel 
                target = "Networks", position = "after") # End InsertTab
      
      i <- i - 1  
    }
    
    removeTab(inputId = "Mainpanels", target = "Data Upload")
    removeTab(inputId = "Network_Tabs", target = "Networks")
    
    
    observeEvent(input$button_treat, {
      
      # Render Download-Buttons
      output$down_smooth_ui <- renderUI({
        div(downloadButton(outputId = "down_smooth", 
                           label = "Download Pre-Processed Data",
                           style = "width:100%;"))
      })
      
      output$down_rep_comb_ui <- renderUI({
        div(downloadButton(outputId = "down_rep_comb", 
                           label = "Download Pooled Replicates",
                           style = "width:100%;"))
      })
      
      # Actual Pre-Processing: Transform data.frame to 4dim-array, applies peak-filtering, smoothing, and normalisation
      # For DataSet A
      my_data_a[is.na(my_data_a)] <- 0
      
      my_data_array_a <- data.array.shiny(x = my_data_a, 
                                        filterpeak = input$box_filterpeak,
                                        normalization = input$box_max, 
                                        lower_bound = input$num_lower, 
                                        smoothing = input$box_loess,
                                        span_value = input$num_span,
                                        names_treatments = metadata_a[[1]],
                                        nr_treatments = metadata_a[[2]],
                                        names_rep = metadata_a[[3]],
                                        nr_rep = metadata_a[[4]],
                                        names_column = metadata_a[[5]],
                                        nr_column = metadata_a[[6]])
      
      my_data_array_a <- my_data_array_a[(rowSums(my_data_array_a) > 0),,, ,drop = FALSE]
      
      # For DataSet B
      my_data_b[is.na(my_data_b)] <- 0
      
      my_data_array_b <- data.array.shiny(x = my_data_b, 
                                          filterpeak = input$box_filterpeak,
                                          normalization = input$box_max, 
                                          lower_bound = input$num_lower, 
                                          smoothing = input$box_loess,
                                          span_value = input$num_span,
                                          names_treatments = metadata_b[[1]],
                                          nr_treatments = metadata_b[[2]],
                                          names_rep = metadata_b[[3]],
                                          nr_rep = metadata_b[[4]],
                                          names_column = metadata_b[[5]],
                                          nr_column = metadata_b[[6]])
      
      my_data_array_b <- my_data_array_b[(rowSums(my_data_array_b) > 0),,, ,drop = FALSE]
    
      # Render Row-Selector UIs
      output$selector_a_ui <- renderUI({
        div(selectInput(inputId = "select_rows_a", 
                        label = "", 
                        choices = rownames(my_data_array_a), 
                        multiple = FALSE, 
                        width = "600px"))
      })
      
      output$selector_b_ui <- renderUI({
        div(selectInput(inputId = "select_rows_b", 
                        label = "", 
                        choices = rownames(my_data_array_b), 
                        multiple = FALSE, 
                        width = "600px"))
      })
      
      output$select_previous_ui <- renderUI({
        actionButton("select_previous", label = HTML("&#9650;"), width = "50px")
        
      })
      
      output$select_next_ui <- renderUI({
        actionButton("select_next", label = HTML("&#9660;"), width = "50px")
        
      })
      
      
#      observeEvent(input$dataset_selector, {
       
        observeEvent(input$select_previous, {
          
          if(input$dataset_selector == "Dataset A"){
            current_a <- which(rownames(my_data_array_a) == input$select_rows_a)
            if(current_a > 1){
              updateSelectInput(session, "select_rows_a",
                                selected = rownames(my_data_array_a)[current_a - 1])
          }
          
          } else if(input$dataset_selector == "Dataset B"){
            current_b <- which(rownames(my_data_array_b) == input$select_rows_b)
            if(current_b > 1){
              updateSelectInput(session, "select_rows_b",
                                selected = rownames(my_data_array_b)[current_b - 1])
          }
          }
        })
      
        observeEvent(input$select_next, {
        
          if(input$dataset_selector == "Dataset A"){
            current_a <- which(rownames(my_data_array_a) == input$select_rows_a)
            if(current_a < length(rownames(my_data_array_a))){
              updateSelectInput(session, "select_rows_a",
                                selected = rownames(my_data_array_a)[current_a + 1])
            }
        
            } else if(input$dataset_selector == "Dataset B"){
              
            current_b <- which(rownames(my_data_array_b) == input$select_rows_b)
            if(current_b < length(rownames(my_data_array_b))){
              updateSelectInput(session, "select_rows_b",
                                selected = rownames(my_data_array_b)[current_b + 1])
            }
            }
        })
        
  #    })
      
      output$down_smooth <- downloadHandler(filename = paste0(Sys.Date(),"_", code,"_PROMISed_1_preprocessed.zip"), 
                                            content = function(file){

                                              data_smooothed_table_a <- data.unarray(x = my_data_array_a, 
                                                                                   names_treatments = metadata_a[[1]], 
                                                                                   names_rep = metadata_a[[3]], 
                                                                                   names_column = metadata_a[[5]])
                                              
                                              data_smooothed_table_b <- data.unarray(x = my_data_array_b, 
                                                                                     names_treatments = metadata_b[[1]], 
                                                                                     names_rep = metadata_b[[3]], 
                                                                                     names_column = metadata_b[[5]])
                                              
                                              pdf(paste0(promised_tmp_wd,"/",
                                                         Sys.Date(),"_", code, "_1_",
                                                         metadata_a[[7]],
                                                         "_RSP_", input$box_filterpeak,
                                                         "_norm_", input$box_max, "_",  input$num_lower*100,
                                                         "_smoothed_", input$box_loess, "_", input$num_span*100,"_", code,
                                                         ".pdf"))
                                              PlotMyArray.shiny(my_data_array_a, names_treatments = metadata_a[[1]])
                                              dev.off()  
                                              
                                              pdf(paste0(promised_tmp_wd,"/",
                                                         Sys.Date(),"_", code,"_1_",
                                                         metadata_b[[7]],
                                                         "_RSP_", input$box_filterpeak,
                                                         "_norm_", input$box_max, "_",  input$num_lower*100,
                                                         "_smoothed_", input$box_loess, "_", input$num_span*100,"_", code,
                                                         ".pdf"))
                                              PlotMyArray.shiny(my_data_array_b, names_treatments = metadata_b[[1]])
                                              dev.off()  
                                              
                                              write.table(data_smooothed_table_a, paste0(promised_tmp_wd,"/",
                                                                                         Sys.Date(),"_", code,"_1_",
                                                                                         metadata_a[[7]],
                                                                                         "_RSP_", input$box_filterpeak,
                                                                                         "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                                         "_smoothed_", input$box_loess, "_", input$num_span*100,"_", code,
                                                                                         ".txt"), sep = "\t")
                                              
                                              write.table(data_smooothed_table_b, paste0(promised_tmp_wd,"/",
                                                                                         Sys.Date(),"_", code,"_1_",
                                                                                         metadata_b[[7]],
                                                                                         "_RSP_", input$box_filterpeak,
                                                                                         "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                                         "_smoothed_", input$box_loess, "_", input$num_span*100,"_", code,
                                                                                         ".txt"), sep = "\t")
                                              
                                              Zip_Files <- list.files(path = promised_tmp_wd,  pattern = code)
                                              zip::zipr(zipfile = file, files = paste0(promised_tmp_wd,"/",Zip_Files))
                                              
                                              file.remove(paste0(path = promised_tmp_wd, "/", list.files(path = promised_tmp_wd, pattern = code)))
                                              
                                            }, contentType = "application/zip")
      
      # 1st Page Mainpanel: 
      # Render Front-Page Plot
      
      observeEvent({input$dataset_selector
                    input$select_rows_a
                    input$select_rows_b}, { # Start observing dataset_selector
      
        if(as.character(input$dataset_selector) == "Dataset A"){
          
          names_treatments  <- metadata_a[[1]]
          nr_treatments     <- metadata_a[[2]]
          names_rep         <- metadata_a[[3]]
          nr_rep            <- metadata_a[[4]]
          names_column      <- metadata_a[[5]]
          my_file_name      <- metadata_a[[7]]
          
          selected_data_array <- my_data_array_a
          row_selector <- input$select_rows_a
          
        } else if(as.character(input$dataset_selector) == "Dataset B"){
          
          names_treatments  <- metadata_b[[1]]
          nr_treatments     <- metadata_b[[2]]
          names_rep         <- metadata_b[[3]]
          nr_rep            <- metadata_b[[4]]
          names_column      <- metadata_a[[5]]
          my_file_name      <- metadata_b[[7]]
          
          selected_data_array <- my_data_array_b
          row_selector <- input$select_rows_b
        }
        
        # Render Checkboxes for Front-Page Plot selection
        output$line_plot_check_ui <- renderUI({
        
          div(checkboxGroupInput(inputId = "line_plot_check", 
                                 label = "Replicates to Plot", 
                                 choiceNames = as.vector(apply(X = as.matrix(names_treatments), 
                                                               MARGIN = 1, 
                                                               FUN = paste, 
                                                               c(1 : as.numeric(nr_rep)), 
                                                               sep = paste0("_", as.character(names_rep), "_"))), 
                                 choiceValues = c(1:c(as.numeric(nr_treatments)*as.numeric(nr_rep))),
                                 selected = c(1:c(as.numeric(nr_treatments)*as.numeric(nr_rep))) ))
        })
      
      
        line_plot_check <- rep(FALSE, times = as.numeric(nr_treatments)*as.numeric(nr_rep))
      
        observeEvent(input$line_plot_check, {
          
          line_plot_check[as.numeric(input$line_plot_check)] <- TRUE
          line_plot_check <- matrix(nrow = as.numeric(nr_treatments),
                                    ncol = as.numeric(nr_rep),
                                    line_plot_check,
                                    byrow = TRUE)
          
          output$plot <- renderPlot({ PlotMyArray.row(x = selected_data_array, 
                                                      my_row_name = row_selector, 
                                                      names_treatments = names_treatments, 
                                                      titel_name = input$select_rows,
                                                      plot_what = line_plot_check)
        })
      })
        
        # End of Page 1
        
        # Start Page 6: Dis-Elution-Score
       
        observeEvent(input$Manhattan_Tabs, {  
          j <- as.numeric(input$Manhattan_Tabs)
          
          # Start reactive manhattan-distances:
          # Calculate manhattan-distances for one selected row, used for boxplots:
          my_manhattan_distances <- reactive({
            
            if(nr_treatments >= 2 & nr_rep >= 3){
              
            anova_input_list <- manhattan.row(x = selected_data_array, 
                                              names_treatments = names_treatments, 
                                              select_rows = row_selector, 
                                              j = as.numeric(input$Manhattan_Tabs) )
            }
          })
         
          # Creating Boxplots   
          
          output[[paste0("boxplot", j)]] <- renderPlot({
            
           if(nr_treatments >= 2 & nr_rep >= 3){
            manhattan_distances <- my_manhattan_distances()
            
            manhattan.anova.boxplot(manhattan_distances = manhattan_distances, 
                                    name_panel = name_panel, 
                                    names_treatments = names_treatments,
                                    j = as.numeric(input$Manhattan_Tabs))
            } else {
        
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 1, paste("Requirements not met: At least two conditions and three replicates required."), 
                   cex = 1.6, col = "black")
              
            }
          })
          
          }) # Stop Observing Manhattan_Tabs
        }) # Stop observing dataset_selector (Pre-Processing)
      
        
          # "Global" Dis-Elution-Score:
          ### Render DES download UI-elements
          
          output$down_table <- downloadHandler(filename = paste0(Sys.Date(),"_", code, "_PROMISed_6_DES.zip"), content = function(file){
             
            # Calcuate DES for both datasets, where possbiel:
            if(metadata_a[[4]] >= 3){
              my_manhattan_a <- manhattan.anova.shiny(my_data_array_a, metadata_a[[1]], pvalue = input$pvalue)
              
              
            } else if(metadata_a[[4]] < 3){
              my_manhattan_a <- matrix(nrow = nrow(my_data_array_a), ncol = 1, "Requirements not met")
              colnames(my_manhattan_a) <- "Dis-Elution-Score"
              
            }
            
            if(metadata_b[[4]] >= 3){
              my_manhattan_b <- manhattan.anova.shiny(my_data_array_b, metadata_b[[1]], pvalue = input$pvalue)
              
            } else if(metadata_b[[4]] < 3){
              my_manhattan_b <- matrix(nrow = nrow(my_data_array_b), ncol = 1, "Requirements not met")
              colnames(my_manhattan_b) <- "Dis-Elution-Score"
              
            }
            
            my_data_table_a <- data.unarray(my_data_array_a, 
                                            metadata_a[[1]], 
                                            metadata_a[[4]], 
                                            metadata_a[[5]])
            colnames(my_data_table_a) <- colnames(my_data_a)
            my_data_table_b <- data.unarray(my_data_array_b, 
                                            metadata_b[[1]], 
                                            metadata_b[[4]], 
                                            metadata_b[[5]]) 
            colnames(my_data_table_b) <- colnames(my_data_b)
            
            my_data_table_a <- cbind(my_manhattan_a, my_data_table_a)
            my_data_table_b <- cbind(my_manhattan_b, my_data_table_b)
            
            write.table(my_data_table_a, 
                        paste0(promised_tmp_wd,"/",
                               Sys.Date(),"_", code,"_6_",
                               metadata_a[[7]],
                               "_RSP_", input$box_filterpeak,
                               "_norm_", input$box_max, "_",  input$num_lower*100,
                               "_smoothed_", input$box_loess, "_", input$num_span*100,
                               "_differential_", input$pvalue*100,"_", code,
                               ".txt")
                        ,sep = "\t")
            
            write.table(my_data_table_b, 
                        paste0(promised_tmp_wd,"/",
                               Sys.Date(),"_", code,"_6_",
                               metadata_b[[7]],
                               "_RSP_", input$box_filterpeak,
                               "_norm_", input$box_max, "_",  input$num_lower*100,
                               "_smoothed_", input$box_loess, "_", input$num_span*100,
                               "_differential_", input$pvalue*100,"_", code,
                               ".txt")
                        ,sep = "\t")
            
            Zip_Files <- list.files(path = promised_tmp_wd,  pattern = code)
            zip::zipr(zipfile = file, files = paste0(promised_tmp_wd,"/",Zip_Files))
            
            file.remove(paste0(path = promised_tmp_wd, "/", list.files(path = promised_tmp_wd, pattern = code)))
            
          })
        
      # Page 2: Combining Replicates
      # Combining Replicates based on Chosen Measurement of Correlation:
      # Combining all Array-Rows
      
      observeEvent({input$cor_method_replicate
                    input$pcc_thresh_replicate
                    input$box_rep
                    input$var_min_peak
                    input$var_limit_small
                    input$var_limit_large}, {
          
                      if(as.numeric(metadata_a[[4]]) == 1){
                        
                        my_data_rep_comb_a <- my_data_array_a
                        
                      } else if(as.numeric(metadata_a[[4]]) > 1){
                        
                        my_data_rep_comb_a <- combine.reps.shiny(my_data_array_a, 
                                                                 PCC_thresh = as.numeric(input$pcc_thresh_replicate), 
                                                                 cor_method = as.character(input$cor_method_replicate),
                                                                 names_treatments = metadata_a[[1]], 
                                                                 rep_box = input$box_rep)
                      }
                      
                      if(as.numeric(metadata_b[[4]]) == 1){
                        
                        my_data_rep_comb_b <- my_data_array_b
                        
                      } else if(as.numeric(metadata_b[[4]]) > 1){
                        
                        my_data_rep_comb_b <- combine.reps.shiny(my_data_array_b, 
                                                                 PCC_thresh = as.numeric(input$pcc_thresh_replicate), 
                                                                 cor_method = as.character(input$cor_method_replicate),
                                                                 names_treatments = metadata_b[[1]], 
                                                                 rep_box = input$box_rep)
                      }
                        
                      #Download Replicate-Combined Data:
                      output$down_rep_comb <- downloadHandler(filename = paste0(Sys.Date(),"_", code, "_PROMISed_2_replicate_pooled.zip"), 
                                                              content = function(file){
                                                                
                                                                data_rep_comb_table_a <- data.unarray(x = my_data_rep_comb_a, 
                                                                                                    names_treatments = metadata_a[[1]], 
                                                                                                    names_rep = metadata_a[[3]], 
                                                                                                    names_column = metadata_a[[5]])
                                                                
                                                                data_rep_comb_table_b <- data.unarray(x = my_data_rep_comb_b, 
                                                                                                      names_treatments = metadata_b[[1]], 
                                                                                                      names_rep = metadata_b[[3]], 
                                                                                                      names_column = metadata_b[[5]])
                                                                
                                                                pdf(paste0(promised_tmp_wd,"/",
                                                                           Sys.Date(),"_", code,"_2_",
                                                                           metadata_a[[7]],
                                                                           "_RSP_", input$box_filterpeak,
                                                                           "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                           "_smoothed_", input$box_loess, "_", input$num_span*100,
                                                                           "_", input$cor_method_replicate, "_", input$pcc_thresh_replicate*100, "_single_", input$box_rep,"_", code,
                                                                           ".pdf"))
                                                                PlotMyArray.shiny(my_data_rep_comb_a, names_treatments = metadata_a[[1]])
                                                                dev.off()  
                                                                
                                                                pdf(paste0(promised_tmp_wd,"/",
                                                                           Sys.Date(),"_", code,"_2_",
                                                                           metadata_b[[7]],
                                                                           "_RSP_", input$box_filterpeak,
                                                                           "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                           "_smoothed_", input$box_loess, "_", input$num_span*100,
                                                                           "_", input$cor_method_replicate, "_", input$pcc_thresh_replicate*100, "_single_", input$box_rep,"_", code,
                                                                           ".pdf"))
                                                                PlotMyArray.shiny(my_data_rep_comb_b, names_treatments = metadata_b[[1]])
                                                                dev.off() 
                                                                
                                                                write.table(data_rep_comb_table_a, paste0(promised_tmp_wd,"/",
                                                                                                          Sys.Date(),"_", code,"_2_",
                                                                                                          metadata_a[[7]],
                                                                                                        "_RSP_", input$box_filterpeak,
                                                                                                        "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                                                        "_smoothed_", input$box_loess, "_", input$num_span*100,
                                                                                                        "_", input$cor_method_replicate, "_", input$pcc_thresh_replicate*100, "_single_", input$box_rep,"_", code,
                                                                                                        ".txt"), sep = "\t")
                                                                
                                                                write.table(data_rep_comb_table_b, paste0(promised_tmp_wd,"/",
                                                                                                          Sys.Date(),"_", code,"_2_",
                                                                                                          metadata_b[[7]],
                                                                                                        "_RSP_", input$box_filterpeak,
                                                                                                        "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                                                        "_smoothed_", input$box_loess, "_", input$num_span*100,
                                                                                                        "_", input$cor_method_replicate, "_", input$pcc_thresh_replicate*100, "_single_", input$box_rep,"_", code,
                                                                                                        ".txt"), sep = "\t")
                                                                
                                                                Zip_Files <- list.files(path = promised_tmp_wd,  pattern = code)
                                                                zip::zipr(zipfile = file, files = paste0(promised_tmp_wd,"/",Zip_Files))
                                                                
                                                                file.remove(paste0(path = promised_tmp_wd, "/", list.files(path = promised_tmp_wd, pattern = code)))
                                                                
                                                                
                                                              }, contentType = "application/zip")
                      
                      observeEvent({input$dataset_selector
                                    input$select_rows_a
                                    input$select_rows_b}, { # Start observing dataset_selector
                          
                          if(as.character(input$dataset_selector) == "Dataset A"){
                            
                            names_treatments  <- metadata_a[[1]]
                            nr_treatments     <- metadata_a[[2]]
                            nr_rep            <- metadata_a[[4]]
                            nr_column         <- metadata_a[[6]]
                            
                            selected_data_array <- my_data_array_a
                            row_selector <- input$select_rows_a
                             
                            selected_data_rep_comb <- my_data_rep_comb_a
                            
                          } else if(as.character(input$dataset_selector) == "Dataset B"){
                            
                            names_treatments  <- metadata_b[[1]]
                            nr_treatments     <- metadata_b[[2]]
                            nr_rep            <- metadata_b[[4]]
                            nr_column         <- metadata_b[[6]]
                            
                            selected_data_array <- my_data_array_b
                            row_selector <- input$select_rows_b
                            
                            selected_data_rep_comb <- my_data_rep_comb_b
                            
                          }
                          
                          # Combine single, selected Row to plot
                            if(as.numeric(nr_rep) == 1){
                              
                              selected_data_rep_comb_row <- selected_data_array[which(rownames(selected_data_array) == row_selector),,, ,drop = FALSE]
                              
                            } else if(as.numeric(nr_rep) > 1){
                              
                              selected_data_rep_comb_row <- combine.reps.shiny(selected_data_array[which(rownames(selected_data_array) == row_selector),,, ,drop = FALSE], 
                                                                         PCC_thresh = as.numeric(input$pcc_thresh_replicate),
                                                                         cor_method = as.character(input$cor_method_replicate),
                                                                         names_treatments = names_treatments, 
                                                                         rep_box = input$box_rep)
                            } 
                          
                          # Plot combined, selected Row:
                          output$repcomb_plot <- renderPlot({
                            
                            PlotMyArray.row(x = selected_data_rep_comb_row, 
                                            my_row_name = row_selector, 
                                            names_treatments = names_treatments, 
                                            titel_name = row_selector,
                                            plot_what = matrix(nrow = as.numeric(nr_treatments), ncol = 1, TRUE))
                          })
                          
                          # Create Heatmap for all Rows of selected data set
                          
                          output$"heatmap_plot" <- renderPlot({
                            
                            my_rep_comb_pheatmap <- matrix(ncol = 0, nrow = nrow(selected_data_rep_comb))
                            
                            k <- 1
                            while(k <= nr_treatments){
                              my_rep_comb_pheatmap <- cbind(my_rep_comb_pheatmap, selected_data_rep_comb[,,k,])
                              k <- k + 1
                            }
                            
                            colnames(my_rep_comb_pheatmap) <- make.unique(colnames(my_rep_comb_pheatmap))
                            my_rep_comb_pheatmap <- my_rep_comb_pheatmap[order(apply(X = my_rep_comb_pheatmap, MARGIN = 1, FUN = which.max)),]
                            
                            my_col_anno <- c()
                            
                            ann_colors2 = list(
                              Treatment = as.character(matrix(nrow = nr_treatments, ncol = 1))
                            )
                            
                            gaps_col <- c()
                            k <- 1
                            while(k <= nr_treatments){
                              my_col_anno    <- c(my_col_anno, rep(names_treatments[k], nr_column)) 
                              #ann_colors2$Treatment[k]  <- pal_startrek(palette = c("uniform"))(7)[k]
                              ann_colors2$Treatment[k] <- c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")[k]
                              gaps_col <- c(gaps_col, k*nr_column)
                              k <- k + 1
                            }
                            
                            names(ann_colors2$Treatment) <- names_treatments
                            
                            my_col_anno <- data.frame(Treatment = my_col_anno)
                            rownames(my_col_anno) <- colnames(my_rep_comb_pheatmap)
                            
                            pheatmap(my_rep_comb_pheatmap, 
                                     show_rownames = FALSE,
                                     show_colnames = FALSE,
                                     cluster_cols  = FALSE,
                                     cluster_rows  = FALSE,
                                     annotation_col = my_col_anno,
                                     cex = 1,
                                     gaps_col = gaps_col,
                                     annotation_legend = TRUE,
                                     annotation_names_col = FALSE,
                                     annotation_colors = ann_colors2,
                                     color = hcl.colors(50, "viridis"))
                            
                          }) # End of PlotOutput: Pheatmap
                        
                          # Plot Deconvoluted
                          observeEvent({input$var_min_peak
                            input$var_limit_small
                            input$var_limit_large
                            input$nodeconvolution_input}, {
                              
                              # For treatments in separate Tabs:
                              observeEvent(input$Deconvo_Tabs, {  
                                
                              
                                if(input$nodeconvolution_input == FALSE){
                                
                              selected_data_rep_comb_row_j <- selected_data_rep_comb_row[,,as.numeric(input$Deconvo_Tabs),, drop = FALSE]
                            
                              g <- matrix(selected_data_rep_comb_row_j, nrow = 1, ncol = dim(selected_data_array)[2])
                              g <- rbind(g, 1, 1)
                              rownames(g) <- c(rownames(selected_data_rep_comb_row_j), "dummy1", "dummy2")
                              g[is.na(g)] <- 0
                              g <- g[(rowSums(g) > 0),, drop = FALSE]
                            
                              selected_data_deconvolute <- data.frame(deconvolute.new(my_data = g, 
                                                              var_min_peak = input$var_min_peak, 
                                                              var_limit_small = input$var_limit_small, 
                                                              var_limit_large = input$var_limit_large))
                           
                              } else if(input$nodeconvolution_input == TRUE){
                                
                                selected_data_deconvolute <- data.frame(matrix(nrow = 0, ncol = c(1+ metadata_a[[6]])))
                                
                                i <- 1
                                while(i <= metadata_a[[2]]){
                                  
                                  selected_data_deconvolute <- rbind(selected_data_deconvolute, data.frame(Treatment = metadata_a[[1]][i], selected_data_rep_comb_row[,,i,]))
                                  
                                  i <- i + 1
                                }
                                
                              }
                           
                              
                            # Render (interactive) Plot showing the deconvoluted Peaks of a selected profile
                              output[[paste0("decon_plot_", as.numeric(input$Deconvo_Tabs))]] <- renderPlot({
                              
                                PlotDeconvoluted(rep_comb = selected_data_rep_comb_row, 
                                                 deconvoluted = selected_data_deconvolute, 
                                                 treatment = as.numeric(input$Deconvo_Tabs), 
                                                 my_row_name = row_selector)
                              
                            })
                            
                              
                              }) # Stop observing Deconvo_Tabs                 
                              
                            }) # Stop Observing Deconvolution-Inputs !  # ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
                            
                        }) # Stop Observing dataset_selector (Replicate Combination)
                      
                      # End of Page 2: Replicate Combination
                      
                      # 3rd Tab: Deconvolution
                        
                      # Plot Deconvoluted
                      observeEvent({input$var_min_peak
                        input$var_limit_small
                        input$var_limit_large
                        input$nodeconvolution_input
                        }, {
                        
                          if(input$nodeconvolution_input == FALSE){
                            
                            my_data_deconvoluted_a <- deconvolute.array(my_data_rep_comb_a, 
                                                                        names_treatment = metadata_a[[1]], 
                                                                        var_min_peak = input$var_min_peak, 
                                                                        var_limit_small = input$var_limit_small, 
                                                                        var_limit_large = input$var_limit_large,
                                                                        name_col = metadata_a[[5]],
                                                                        nr_col = metadata_a[[6]])                      
                            
                            my_data_deconvoluted_b <- deconvolute.array(my_data_rep_comb_b, 
                                                                        names_treatment = metadata_b[[1]], 
                                                                        var_min_peak = input$var_min_peak, 
                                                                        var_limit_small = input$var_limit_small, 
                                                                        var_limit_large = input$var_limit_large,
                                                                        name_col = metadata_b[[5]],
                                                                        nr_col = metadata_b[[6]])     
                        
                          } else if(input$nodeconvolution_input == TRUE){

                            my_data_deconvoluted_a <- numeric()
                            my_data_deconvoluted_b <- numeric()
                            
                            i <- 1
                            while(i <= metadata_a[[2]]){
                              
                              nondecon_data_a <- (as.matrix(my_data_rep_comb_a[,,i,]))
                              nondecon_data_a <- data.frame(Treatment = metadata_a[[1]][i], nondecon_data_a)
                              
                              my_data_deconvoluted_a <- rbind(my_data_deconvoluted_a, nondecon_data_a)
                                                       
                              nondecon_data_b <- (as.matrix(my_data_rep_comb_b[,,i,]))
                              nondecon_data_b <- data.frame(Treatment = metadata_b[[1]][i], nondecon_data_b)
                              
                              my_data_deconvoluted_b <- rbind(my_data_deconvoluted_b, nondecon_data_b)
                              
                              i <- i + 1
                            }
                           
                            empty_profiles_a <- rowSums(my_data_deconvoluted_a[,-1]) == 0 
                            empty_profiles_b <- rowSums(my_data_deconvoluted_b[,-1]) == 0
                           
                            my_data_deconvoluted_a <- my_data_deconvoluted_a[!empty_profiles_a,]
                            my_data_deconvoluted_b <- my_data_deconvoluted_b[!empty_profiles_b,]
                            
                           }
                        
                        output$down_decon <- downloadHandler(filename = paste0(Sys.Date(),"_", code, "_PROMISed_3_deconvoluted.zip"), 
                                                                content = function(file){
                                                                  
                                                                  write.table(my_data_deconvoluted_a, paste0(promised_tmp_wd,"/",
                                                                                                             Sys.Date(),"_", code,"_3_",
                                                                                                             metadata_a[[7]],
                                                                                                             "_RSP_", input$box_filterpeak,
                                                                                                             "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                                                             "_smoothed_", input$box_loess, "_", input$num_span*100,
                                                                                                             "_", input$cor_method_replicate, "_", input$pcc_thresh_replicate*100, "_single_", input$box_rep,
                                                                                                             "_decon_", input$var_limit_small*100, "_", input$var_min_peak*100, "_", input$var_limit_large*100,"_", code,
                                                                                                             ".txt"), sep = "\t")
                                                                    
                                                                  write.table(my_data_deconvoluted_b, paste0(promised_tmp_wd,"/",
                                                                                                             Sys.Date(),"_", code,"_3_",
                                                                                                             metadata_b[[7]],
                                                                                                             "_RSP_", input$box_filterpeak,
                                                                                                             "_norm_", input$box_max, "_",  input$num_lower*100,
                                                                                                             "_smoothed_", input$box_loess, "_", input$num_span*100,
                                                                                                             "_", input$cor_method_replicate, "_", input$pcc_thresh_replicate*100, "_single_", input$box_rep,
                                                                                                             "_decon_", input$var_limit_small*100, "_", input$var_min_peak*100, "_", input$var_limit_large*100,"_", code,
                                                                                                             ".txt"), sep = "\t")
                                                                  
                                                                  Zip_Files <- list.files(path = promised_tmp_wd,  pattern = code)
                                                                  zip::zipr(zipfile = file, files = paste0(promised_tmp_wd,"/",Zip_Files))
                                                                  
                                                                  file.remove(paste0(path = promised_tmp_wd, "/", list.files(path = promised_tmp_wd, pattern = code)))
                                                                  
                                                                }, contentType = "application/zip")
                        
                        # Start Data-Integration:
                        observeEvent(input$button_integration, {
                          
                          decon_list_a <- deconvoluted.list(my_data_deconvoluted_a)  
                          decon_list_b <- deconvoluted.list(my_data_deconvoluted_b)
                          
                          
                          # Download Correlation Tables:
                          output$down_corr_table_zip <- downloadHandler(filename = paste0(Sys.Date(),"_", code,"_PROMISed_4_Cor_Tables.zip"), 
                                                                        content = function(file){
                            
                            #setwd(promised_tmp_wd)
                            
                            m <- 1
                            while(m <= length(decon_list_a)){
                              
                              x <- as.matrix(decon_list_a[[m]])
                              y <- as.matrix(decon_list_b[[m]])
                              
                              my_correl_table <- cor(x = t(x), y = t(y), method = as.character(input$cor_method_integration))
                              
                              write.table(my_correl_table, paste0(promised_tmp_wd,"/",
                                                                  Sys.Date(),"_", code,"_4_",
                                                                  metadata_a[[7]], 
                                                                  "_", metadata_b[[7]],
                                                                  "_", as.character(input$cor_method_integration), 
                                                                  "_cor_table_", 
                                                                  code,"_", 
                                                                  metadata_a[[1]][m], 
                                                                  ".txt"), sep = "\t")
                              
                              m <- m + 1  
                            }
                            
                            Zip_Files <- list.files(path = promised_tmp_wd,  pattern = code)
                            zip::zipr(zipfile = file, files = paste0(promised_tmp_wd,"/",Zip_Files))
                            
                            file.remove(paste0(path = promised_tmp_wd, "/", list.files(path = promised_tmp_wd, pattern = code)))
                            
                             
                          }, contentType = "application/zip")
                          
                          # Create Plots while observing active Panels:
                          
                          observeEvent({input$Integration_Tabs
                                        input$dataset_selector
                                        input$select_rows_a
                                        input$select_rows_b}, {
                            
                            # Subset Data on active Tab
                            j <- as.numeric(input$Integration_Tabs)
                            my_data_decon_a_j <- data.frame(decon_list_a[[j]])
                            my_data_decon_b_j <- data.frame(decon_list_b[[j]]) 
                            
                            # Create Targeted plot using BOTH datasets and selecotrs
                            output[[paste0("target_plot_", j)]] <- renderPlot({
                              
                              Plot2selections(data1 = my_data_decon_a_j,
                                              data2 = my_data_decon_b_j,
                                              selector1 = input$select_rows_a,
                                              selector2 = input$select_rows_b)
                              
                            })
                            
                            # Create Plots requiring data-selection: Intersection and "Best Hit"
                         
                                # Select Deconvoluted Data
                                if(input$dataset_selector == "Dataset A"){
                                  
                                  selected_decon_data <- my_data_decon_a_j
                                  other_decon_data    <- my_data_decon_b_j
                                  
                                  selected_decon_list <- decon_list_a
                                  other_decon_list    <- decon_list_b
                                  names_treatments    <- metadata_a[[1]]
                                  nr_treatments       <- metadata_a[[2]]
                                  
                                  best_hit_selector   <- input$select_rows_a
                                  
                                } else if(input$dataset_selector == "Dataset B"){
                                    
                                  selected_decon_data <- my_data_decon_b_j
                                  other_decon_data    <- my_data_decon_a_j
                                  
                                  
                                  selected_decon_list <- decon_list_b
                                  other_decon_list    <- decon_list_a
                                  names_treatments    <- metadata_b[[1]]
                                  nr_treatments       <- metadata_b[[2]]
                                   
                                  best_hit_selector   <- input$select_rows_b
                                }   
                                
                                output[[paste0("hits_plot_", j)]] <- renderPlot({
                                  
                                  PlotBestHit(decon_data1 = selected_decon_data, 
                                              decon_data2 = other_decon_data, 
                                              selector = best_hit_selector, 
                                              pcc_threshold = input$pcc_slider_integration,
                                              method = as.character(input$cor_method_integration))
                                })
                                
                                cor_display <- cor_reactive(data1 = selected_decon_data,
                                                            data2 = other_decon_data , 
                                                            pattern = best_hit_selector, 
                                                            method = as.character(input$cor_method_integration))
                                
                                output[[paste0("pcc_table_", j)]]  <- DT::renderDataTable(    
                                  
                                  formatRound(datatable(cor_display), 
                                              columns = c(1:ncol(cor_display)),
                                              digits = 3),
                                              
                                  options = list(
                                    dom = 'Bfrltip',
                                    pageLength = 5,
                                    scrollX = TRUE,
                                    scrollCollapse = TRUE,
                                    searchHighlight = TRUE
                                  ),
                                  style = "default"
                                )
                                
                              if(input$Integration_Tabs == "Intersections of Conditions"){
                                
                                Coelution_list <- ListCoelutions(data1 = selected_decon_list,
                                                                 data2 = other_decon_list,
                                                                 selector = best_hit_selector,
                                                                 pcc_threshold = input$pcc_slider_integration,
                                                                 method = input$cor_method_integration,
                                                                 names_treatments = names_treatments,
                                                                 nr_treatments = nr_treatments)()
                                
                                output$down_overlap <- downloadHandler(filename = paste0("Sections_", 
                                                                                         best_hit_selector,"_", 
                                                                                         as.character(input$cor_method_integration),"_", 
                                                                                         input$pcc_slider_integration, ".txt"), 
                                                                         content = function(file) {
                                                                           
                                                                           intersections_temp <- GetIntersections(x = Coelution_list) 
                                                                           
                                                                           write.table(rbind(colSums(!is.na(intersections_temp)), intersections_temp), file, sep = "\t")
                                                                           
                                                                         })
                                observeEvent(input$euler_input, {
                                output$venn_set <- renderPlot({ 
                                
                                  venn_colors <- list("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                  
                                  if(length(Coelution_list) != 0){
                                  
                                    if(input$euler_input == "Venn Diagram"){
                                  
                                      if(length(Coelution_list) > 1){
                                        
                                        plot(venn(Coelution_list),
                                             fills = list(fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[1:length(Coelution_list)], alpha = 0.7),
                                             quantities = list( cex = 2),
                                             names = list(cex = 3),
                                             legend = list(labels = names(Coelution_list), cex = 2))
                                 
                                      } else {
                                        
                                        Coelution_list2 <- Coelution_list[lapply(Coelution_list, length)>0]
                                        venn_colors <- venn_colors[1:length(Coelution_list)]
                                        venn_colors <- venn_colors[lapply(Coelution_list, length)>0]
                                        
                                        plot(euler(Coelution_list2),
                                             fills = list(fill = unlist(venn_colors), alpha = 0.7),
                                             quantities = list( cex = 2),
                                             names = list(cex = 3),
                                             legend = list(labels = names(Coelution_list2), cex = 2))
                                        
                                      }
                                        
                                    } else if (input$euler_input == "Euler Diagram"){
                                      
                                      Coelution_list2 <- Coelution_list[lapply(Coelution_list, length)>0]
                                      venn_colors <- venn_colors[lapply(Coelution_list, length)>0]
                                      
                                      plot(euler(Coelution_list2),
                                           fills = list(fill = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[1:length(Coelution_list)], alpha = 0.7),
                                           quantities = list( cex = 2),
                                           names = list(cex = 3),
                                           legend = list(labels = names(Coelution_list2), cex = 2))
                                      
                                    }
                                  } else {
                                    
                                    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                                    text(x = 0.5, y = 1, paste("No co-fractionating molecules found using curent setting."), 
                                         cex = 1.6, col = "black")
                                    
                                  }
                                })
                                
                                })

                              } # End of "Intersections of Conditions"   
                                
                              }) # Stop observing Data-Selections # and Integration_Tabs
                          # End of Page 4 Integration  
                          
                          # Start of Page 5 Network Analysis:
                          
                          observeEvent({input$FirstOrder
                                        input$dataset_selector
                                        input$select_rows_a
                                        input$select_rows_b
                                        input$network_filter}, {
                           
                                          observeEvent(input$Network_Tabs, {
                                            
                                            nettab <- as.numeric(as.numeric(input$Network_Tabs))
                                            my_data_decon_a <- data.frame(decon_list_a[[nettab]])
                                            my_data_decon_b <- data.frame(decon_list_b[[nettab]]) 
                                            
                                           
                                          if(input$network_filter == "network_filter_0"){
                                            
                                            my_cor_table3 <- cor(x = t(my_data_decon_b), y = t(my_data_decon_a), method = as.character(input$cor_method_integration))
                                            
                                          } else if (input$network_filter == "network_filter_1"){
                                            
                                            my_data_decon_a_sel <- my_data_decon_a[grep(pattern = input$select_rows_a, rownames(my_data_decon_a), fixed = TRUE),, drop = FALSE]
                                            
                                            if(nrow(my_data_decon_a_sel != 0)){  
                                            
                                              my_cor_table3 <- cor(x = t(my_data_decon_b), y = t(my_data_decon_a_sel), method = as.character(input$cor_method_integration))
                                            
                                             } else {
                                              
                                              my_cor_table3 <- data.frame(EntryNotFound = c(0,1),
                                                                          TrySomethingElse = c(1,0))
                                              
                                              colnames(my_cor_table3) <- rownames(my_cor_table3)
                                              
                                             }
                                            
                                          } else if (input$network_filter == "network_filter_2"){
                                            
                                            my_data_decon_b_sel <- my_data_decon_b[grep(pattern = input$select_rows_b, rownames(my_data_decon_b), fixed = TRUE),, drop = FALSE]
                                            
                                            if(nrow(my_data_decon_b_sel != 0)){  
                                              
                                              my_cor_table3 <- cor(x = t(my_data_decon_b_sel), y = t(my_data_decon_a), method = as.character(input$cor_method_integration))
                                          
                                            } else {
                                              
                                              my_cor_table3 <- data.frame(EntryNotFound = c(0,1),
                                                                          TrySomethingElse = c(1,0))
                                            
                                              colnames(my_cor_table3) <- rownames(my_cor_table3)
                                                
                                            }
                                          }
                            
                                            my_cor_table3[my_cor_table3 < as.numeric(input$pcc_slider_integration)] <- 0
                                            diag(my_cor_table3) <- 0
                                            
                                            ##Networking
                                            ### IGRAPH NETWORK ###
                                            
                                            network_raw <- custom_igraph(my_cor_table3, pcc = as.numeric(input$pcc_slider_integration))
                                            network_igraph <- network_raw
                                            network_decomposed <- decompose(network_raw, "strong")
                                            
                                            if(length(network_decomposed) > 1){
                                              network_igraph <- network_decomposed[[1]]
                                              for(i in c(2:length(network_decomposed))){
                                                network_igraph <- network_igraph + network_decomposed[[i]]
                                              }
                                            } else {
                                              network_igraph <- network_decomposed[[1]]
                                            }
                                            
                                            my_clust <- cluster_louvain(network_igraph, weights = E(network_igraph)$weight)
                                            
                                            coreness <- coreness(network_igraph)
                                            
                                            coreness_scale <- scale_coreness(coreness)[[1]]
                                            coreness_color <- scale_coreness(coreness)[[2]]
                                            
                                            output[[paste0("network_files_", nettab)]] <- downloadHandler(filename = paste0(Sys.Date(),"_", code,"_PROMISed_5_Network_tables_", metadata_a[[1]][nettab],
                                                                                                                            "_pcc_", input$pcc_slider_integration, 
                                                                                                                            ".zip"),
                                                                                                          content = function(file) {

                                                                                                            write.table(igraph::as_data_frame(network_igraph, what = "edges"),
                                                                                                                        paste0(promised_tmp_wd,"/",
                                                                                                                               Sys.Date(),"_", code,
                                                                                                                               "_5_",metadata_a[[7]], "_" ,
                                                                                                                               metadata_b[[7]], "_", metadata_a[[1]][nettab], 
                                                                                                                               "_", as.character(input$cor_method_integration), 
                                                                                                                               "_", input$pcc_slider_integration,
                                                                                                                               "_network_edges",".txt"), 
                                                                                                                        sep = "\t")
                                                                                                            write.table(cbind(igraph::as_data_frame(network_igraph, what = "vertices"), my_clust$membership, coreness), 
                                                                                                                        paste0(promised_tmp_wd,"/",
                                                                                                                               Sys.Date(),"_", code,
                                                                                                                               "_5_",metadata_a[[7]], "_" ,
                                                                                                                               metadata_b[[7]], "_", metadata_a[[1]][nettab], 
                                                                                                                               "_", as.character(input$cor_method_integration), 
                                                                                                                               "_", input$pcc_slider_integration,
                                                                                                                               "_network_nodes",".txt"), 
                                                                                                                        sep = "\t")
                                                                                                            
                                                                                                            Zip_Files <- list.files(path = promised_tmp_wd,  pattern = code)
                                                                                                            zip::zipr(zipfile = file, files = paste0(promised_tmp_wd,"/",Zip_Files))
                                                                                                            
                                                                                                            file.remove(paste0(path = promised_tmp_wd, "/", list.files(path = promised_tmp_wd, pattern = code)))
                                                                                                            
                                                                                                          })
                                            
                                            ##### VISNETWORK #####
                                            
                                            recomposed_layout <- list(layout_(network_igraph, with_fr(), component_wise()),
                                                                      layout_(network_igraph, in_circle(), component_wise()),
                                                                      layout_(network_igraph, as_tree(), component_wise()),
                                                                      layout_(network_igraph, randomly(), component_wise()))
                                              
                                            observeEvent(input$node_color_input, {
                                              
                                              mycol_network <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                                 "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                                                 "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                                                                 "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                                                 "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                                                 "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
                                              
                                              legendary_node_list <- list(data.frame(label = sort(unique(my_clust$membership)),
                                                                                     color = mycol_network[sort(unique(my_clust$membership))]),
                                                                          data.frame(label = format(seq(from = min(coreness), to = max(coreness), length.out = length(unique(coreness_scale))), digits = 1),
                                                                                     color = coreness_color[sort(unique(coreness_scale))]))
                                              
                                              colouring_list <- list(mycol_network[my_clust$membership], coreness_color[coreness_scale])
                                              
                                              V(network_igraph)$color <- colouring_list[[as.numeric(input$node_color_input)]] 
                                              V(network_igraph)$cluster <- my_clust$membership
                                              V(network_igraph)$coreness <- coreness_scale
                                              E(network_igraph)$color <- "black"
                                              
                                              network_vis <- toVisNetworkData(network_igraph)
                                              
                                              observeEvent(input$my_layout, {
                                                
                                                output[[paste0("visnetwork_network_", nettab)]] <- renderVisNetwork({  #   #FeatureStacking
                                                  
                                                  visIgraph(network_igraph) %>%
                                                    
                                                    visIgraphLayout(layout  = "layout.norm",
                                                                    physics = FALSE,
                                                                    smooth  = FALSE,
                                                                    layoutMatrix = recomposed_layout[[as.numeric(input$my_layout)]]
                                                    ) %>%
                                                    
                                                    visPhysics(stabilization = FALSE) %>%
                                                    
                                                    visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
                                                               manipulation = FALSE, 
                                                               selectedBy = "cluster",
                                                               nodesIdSelection = list(enabled = TRUE)  
                                                    ) %>%  
                                                    
                                                    visInteraction(navigationButtons = TRUE,
                                                                   keyboard = TRUE,
                                                                   multiselect = TRUE,
                                                                   selectConnectedEdges = TRUE,
                                                                   dragNodes = TRUE,
                                                                   zoomView = TRUE
                                                    ) %>%
                                                    
                                                    visLegend(useGroups = FALSE,
                                                              addNodes = legendary_node_list[[as.numeric(input$node_color_input)]],
                                                              main = c("Community", "k-coreness")[[as.numeric(input$node_color_input)]],
                                                              position = "right", ncol = 2, zoom = FALSE) %>%
                                                    
                                                    visEdges(color = "black", width = 1) %>%
                                                    visNodes(borderWidth = 1)
                                                  
                                                })
                                              })
                                              
                                            })    
                                          }) # Stop observing Network Tabs
                            
                          }) # Stop observing Network Input
                          
                        }) # Stop observing button_integration
                        }) # Stop Observing Deconvolution-Inputs !  # ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
                        }) # Stop observing Rep-Comb Inputs
        
  }) # Stop observing button_treat
  }) # Stop observing upload_button
  }) # Stop observing Uploaded Files A and B

  # Load Title-Image
  output$my_image_processing <- renderImage({
    
    list(src = "./www/PROMISed APP Logo.png",
         contentType = "image/png",
         width = 2309/5,
         height = 529/5,
         alt = "Could not find PROMISed Logo"
    )
    
  }, deleteFile = FALSE)
  
  output$logo_mpimp <- renderImage({
    
    list(src = "./www/mpimptext_english_green_transparent_0.png",
         contentType = "image/png",
         width = 2362/10,
         height = 1116/10,
         alt = "Could not find MPIMP Golm Logo"
    )
    
  }, deleteFile = FALSE)
  
  output$logo_minerva <- renderImage({
    
    list(src = "./www/logo_mpg_minerva.png",
         contentType = "image/png",
         width = 2362/20,
         height = 2362/20,
         alt = "Could not find MPIMP Golm Logo"
    )
    
  }, deleteFile = FALSE)
  
  
  
  # Provide Demo Data to download:
  output$download_demo_data_ui <- renderUI({
    div(
      h4("Download Demo-Datasets: "),
      downloadButton(outputId = "download_demo_metabolites", 
                     label = "Metabolites",  
                     style = "width:300px;"),
      br(""),
      downloadButton(outputId = "download_demo_proteins", 
                     label = "Proteins",  
                     style = "width:300px;")
    )
  })
  
  output$download_demo_metabolites <- downloadHandler(filename = "Metabolite_Demodata_PROMISed.txt", content = function(file){
    write.table(metabolite_demo_data, file, sep = "\t")
  })
  
  output$download_demo_proteins <- downloadHandler(filename = "Protein_Demodata_PROMISed.txt", content = function(file){
    write.table(protein_demo_data, file, sep = "\t")
  })
  
  
} # End of SERVER Mainframe
  
shinyApp(ui = ui, server = server)
