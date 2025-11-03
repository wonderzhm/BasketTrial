library(shiny)
library(rhandsontable)

shinyUI(
  navbarPage(title = tags$div(
    tags$img(src="https://i.ibb.co/yFzV7v9r/IDSWG-logo.png", height="40px"), 
    "Basket Trial Design"),
             
             # ---- Home tab ----
             tabPanel("Home",
                      tags$div(
                        img(src = "https://i.ibb.co/60Xp19rS/bskettrial.png", width = "30%"),
                        style = "text-align: center; margin: 0 auto;"
                      ),
                      h2("Basket Trial Design Using Local Power Prior"),
                      p("In oncology drug development, the traditional approach is to evaluate a novel drug's effectiveness one tumor at a time, which is an inefficient process and causes delays in
lack of landscape view of a drug's profile in cancers. Bayesian basket designs, on the other hand, are innovative approaches that allow simultaneous evaluation of a novel
drug across multiple tumor types."),
                      p("The recently proposed local power prior approach is a Bayesian method that enhances the statistical inference by dynamically borrowing information across various
cohorts, while effectively controls the bias. Although computationally convenient, users still face challenges how to use the published R package by Zhou et al (2025)."),
                      p("In this project, we developed an interactive R Shiny app that implements the proposed local power prior method, JSD method, and no-borrowing method (i.e. independent
cohorts) with a graphical interface. This app allows users to specify trial parameters such as sample size, number of cohorts, and null or alternative hypotheses based on
response rate. The operating characteristics such as basket-wise efficacy cutoffs, posterior probabilities, and success rates, can be explored by simulations under various
design options. So statistician practitioners can optimize their study design very conveniently."),
                      p("In addition, the R shiny app can also perform analysis for a basket clinical trial. The statistical analysis results include:"),
                      tags$ul(
                        tags$li("the posterior distribution of response rate for each tumor type and graphically displayed."),
                        tags$li("the posterior mean, median, and credible intervals for each tumor type."),
                        tags$li("the statistical conclusion whether a tumor type considered promising."),
                      ),
             p("Reference: "), 
             tags$ul(
               tags$li("Zhou H, Shen R, Wu S and He P. (2025). A Bayesian Basket Trial Design Using Local Power Prior. Biometrical Journal. DOI: 10.1002/bimj.70069.")
             ),
             p("App Developer and Maintainer:"),
             tags$ul(
               tags$li("Sai Suman Ravi, ss4555-AT-scarletmail.rutgers.edu")
             ),
             p("***The App was created under IDSWG Oncology WG (https://oncologytrialdesign.org/). MIT License***")
             ),
                 # ---- Trial Design tab ----
             tabPanel("Trial Design",
                      sidebarLayout(
                        sidebarPanel(
                          numericInput("p0", "Null Response Rate (p0)", value = 0.15, min = 0, max = 1, step = 0.01),
                          numericInput("p1", "Target Response Rate (p1)", value = 0.30, min = 0, max = 1, step = 0.01),
                          numericInput("alpha", "Basket-wise Type I Error Rate (one-sided)", value = 0.10, min = 0, max = 1, step = 0.01),
                          
                          numericInput("num_baskets", "Number of Baskets (B)", value = 5, min = 1),
                          uiOutput("basket_sample_sizes"),
                          
                          numericInput("num_interim", "Number of Interim Futility Analyses (K)", value = 0, min = 0),
                          uiOutput("interim_inputs"),
                          
                          numericInput("ntrial_design", "Number of Simulations for Cutoffs", value = 5000, min = 100, max = 50000),
                          
                          selectInput("design_method", "Design Method",
                                      choices = c("Independent", "local-PP-PEB", "local-PP-GEB", "JSD"),
                                      selected = "Independent"
                          ),
                          
                          # Beta prior parameters (always shown)
                          numericInput("b1", "Beta Prior a0", value = 0.15, min = 0),
                          numericInput("b2", "Beta Prior b0", value = 0.85, min = 0),
                          
                          # Method-specific parameters
                          conditionalPanel(
                            condition = "input.design_method == 'local-PP-PEB' || input.design_method == 'local-PP-GEB'",
                            numericInput("a", "Borrowing Parameter (a)", value = 1, min = 0),
                            numericInput("delta", "Similarity Threshold (delta)", value = 0.4, min = 0, max = 1, step = 0.01)
                          ),
                          
                          conditionalPanel(
                            condition = "input.design_method == 'JSD'",
                            numericInput("epsilon", "Epsilon Parameter", value = 2, min = 0),
                            numericInput("tau", "Tau Threshold", value = 0.3, min = 0, max = 1, step = 0.01)
                          ),
                          
                          actionButton("calculate_cutoffs", "Calculate Basket-wise Efficacy Cutoffs", 
                                       class = "btn-primary")
                        ),
                        
                        mainPanel(
                          h3("Basket-wise Efficacy Cutoffs"),
                          tableOutput("efficacy_cutoff_table")
                        )
                      )
             ),
             
             # ---- Simulation tab ----
             tabPanel("Simulation",
                      sidebarLayout(
                        sidebarPanel(
                          numericInput("num_simulations", "Number of Simulations", value = 1000, min = 1, max = 50000),
                          numericInput("set_seed", "Set Seed (optional)", value = 1024),
                          
                          uiOutput("simulation_scenarios_ui"),
                          
                          br(),
                          actionButton("run_simulation", "Run Simulation", class = "btn-success")
                        ),
                        
                        mainPanel(
                          h3("Simulation Results"),
                          uiOutput("simulation_progress"),
                          tableOutput("simulation_results_table1"),
                          uiOutput("simulation_results_table2")
                        )
                      )
             ),
             
             # ---- Trial Analysis tab ----
             tabPanel("Trial Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Trial Data Input"),
                          numericInput("num_baskets_analysis", "Number of Baskets", value = 5, min = 1, max = 20),
                          
                          hr(),
                          h5("Enter observed data for each basket:"),
                          uiOutput("trial_analysis_table_inputs"),
                          
                          br(),
                          actionButton("analyze_trial", "Analyze the Trial", class = "btn-primary")
                        ),
                        
                        mainPanel(
                          h3("Trial Analysis Results"),
                          p("This analysis uses the design method and parameters from the Trial Design tab."),
                          tableOutput("trial_analysis_results")
                        )
                      )
             )
  )
)