library(shiny)
library(cowplot) # for plot_grid
library(ggplot2)
library(ggpubr)
library("effsize") # inc cohen.d, wasn't working with "psych" package on shinyapps.io
library("MASS")
library(pwr)
library(lme4)

# user define function to run simulation ----
simulate <- function(input){
  
  # get input
  input_gp1_mean <- input$input_mean
  input_gp2_mean <- 0 ############## could also make this an input
  input_gp1_sd <- input$input_sd # true SD
  input_gp2_sd <- input$input_sd ############## could also make this different
  N <- input$N
  input_noise_sd <- input$input_noise_sd
  
  # run calculations
  true_gp1 <- rnorm(n=N, mean=input_gp1_mean, sd=input_gp1_sd)
  true_gp2 <- rnorm(n=N, mean=input_gp2_mean, sd=input_gp2_sd)
  
  obs_gp1_t1 <- true_gp1 + rnorm(n=N, mean=0, sd=input_noise_sd)
  obs_gp2_t1 <- true_gp2 + rnorm(n=N, mean=0, sd=input_noise_sd)
  obs_gp1_t2 <- true_gp1 + rnorm(n=N, mean=0, sd=input_noise_sd)
  obs_gp2_t2 <- true_gp2 + rnorm(n=N, mean=0, sd=input_noise_sd)        
  
  sim.df <- data.frame( true = c(true_gp1, true_gp2), 
                        obs_t1 = c(obs_gp1_t1, obs_gp2_t1),
                        obs_t2 = c(obs_gp1_t2, obs_gp2_t2),
                        group = c(rep("Group 1",N), rep("Group 2",N))
  ) 
  
  true_pval <- t.test(sim.df[which(sim.df$group=="Group 1"), "true"], sim.df[which(sim.df$group=="Group 2"), "true"])$p.value
  t1_pval <- t.test(sim.df[which(sim.df$group=="Group 1"), "obs_t1"], sim.df[which(sim.df$group=="Group 2"), "obs_t1"])$p.value
  t2_pval <- t.test(sim.df[which(sim.df$group=="Group 1"), "obs_t2"], sim.df[which(sim.df$group=="Group 2"), "obs_t2"])$p.value
  
  # observed t1/t2 correlation (or ICC)
  c1 <- cor(sim.df[which(sim.df$group=="Group 1"), "obs_t1"], sim.df[which(sim.df$group=="Group 1"), "obs_t2"])
  c2 <- cor(sim.df[which(sim.df$group=="Group 2"), "obs_t1"], sim.df[which(sim.df$group=="Group 2"), "obs_t2"])
  # # ICC(cbind(obs_gp2_t1, obs_gp2_t2)) # or could use this
  
  # effect sizes
  d_true <- cohen.d(sim.df$true, sim.df$group)$estimate
  d_t1 <- cohen.d(sim.df$obs_t1, sim.df$group)$estimate
  d_t2 <- cohen.d(sim.df$obs_t2, sim.df$group)$estimate
  
  # correlation between actual and sampled value
  acc1 <- cor(sim.df$obs_t1, sim.df$true)
  acc2 <- cor(sim.df$obs_t2, sim.df$true)
  
  stats.df <- data.frame( true_pval=true_pval, t1_pval=t1_pval, t2_pval=t2_pval, c1=c1, c2=c2, d_true=d_true, d_t1=d_t1, d_t2=d_t2, acc1=acc1, acc2=acc2 )  
  
  return(list(sim=sim.df, stats=stats.df))
  
}

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(h1("Reliability and Group Differences Simulation", align = "center"), windowTitle="Group Differences Simulation"),
  
  fluidRow(
    column(10, offset=1,
           p(style = "text-align:center;", "With this app, you can simulate an experiment to test for group differences, with the goal of illustrating how measurement error leading to low individual reliability can influence power to detect group differences.  The simulation allows you to explore how factors like expected magnitude of the group difference, population variation, sample size, measurement error, and reliability interact."),
           p(style = "text-align:center;", "Results are given for a single simulation to allow for visualization of an example data set, as well as distributions over 100 observations to illustrate central tendencies for a given parameter set."),
           p(style = "text-align:center;", "Experiment with the parameters on the left, and view the results in the middle.  Press the \"Re-run simulation\" button to simulate new experiments with the same parameters, and watch how the results can vary.  Or, press one of the buttons on the right corresponding to a common scenario to update the parameters accordingly.")
    )
  ),
  
  hr(),
  
  # Sidebar layout with input and output definitions ----
  fluidRow(
    
    # Sidebar panel for inputs ----
    column(3, style = "background-color:#f4f4f4;",
           
           # Input sliders
           h4("Parameters"),
           sliderInput(inputId = "input_mean", label = "Expected Group Difference", min = 0, max = 10, value = 1),
           sliderInput(inputId = "input_sd", label = "True Population Variation", min = 0.01, max = 10, value = 3),
           sliderInput(inputId = "N", label = "Sample Size", min = 10, max = 1500, value = 50),
           sliderInput(inputId = "input_noise_sd", label = "Measurement Error", min = .01, max = 10, value = 3),
           actionButton("run", label = "Re-run simulation"),
           htmlOutput("text_params"),
           br()
           
    ),
    
    # Main panel for displaying outputs ----
    column(7, 
           
           h4("Results - single simulation of two time points"),
           br(),
           plotOutput(outputId = "distPlot"),
           br(),
           h4("Results - distributions over 100 simulations of two time points"),
           br(),
           plotOutput(outputId = "distPlot100")
           
    ),
    
    # Sidebar panel (R) ----
    column(2, style = "background-color:#f4f4f4;",
           
           h4("Example scenarios"),
           p(strong("1: Robust, Unreliable, Underpowered."),  "A common scenario where groups truly are different, but measurement error leads to low individual reliability, limited ability to detect group differences, and underestimates of effect sizes.  Use the slider to decrease measurement error and notice how this increases reliability as well as the number of times group differences are detected."),
           actionButton("RobUnrelUnpow", label = "Run this scenario"),
           br(),
           br()
           # p(strong("2: Reliable, Well-powered."),  "With a reliable measure and lots of subjects, you can be well-powered to detect associations with another variable even if you don't observe group-level activation."),
           # actionButton("RelPow", label = "Run this scenario"),
           # br(),
           # br()
           
    )
    
  )
)


# Define server logic required to create outputs ----
server <- function(input, output, session) {
  
  vals <- reactiveValues()
  observe({
    
    input$run # access re-run button, so that if it's clicked this will run again
    
    # run 1 iteration
    result <- simulate(input)
    
    # also run 100 iterations
    stats100.df <- data.frame( )
    for( i in 1:100 ){
      stats.df <- simulate(input)$stats
      stats100.df <- rbind(stats100.df, stats.df)
    }
    
    # save variables
    vals$sim.df <- result$sim 
    vals$stats.df <- result$stats     
    vals$stats100.df <- stats100.df
    
  })
  
  # Update for different example scenarios ----
  observeEvent(input$RobUnrelUnpow, {
    updateSliderInput(session, "input_mean", value = 3)
    updateSliderInput(session, "input_sd", value = 3)
    updateSliderInput(session, "N", value = 50)
    updateSliderInput(session, "input_noise_sd", value = 9)
  })
  
  # observeEvent(input$RelPow, {
  #   updateSliderInput(session, "input_mean", value = 0)
  #   updateSliderInput(session, "input_sd", value = 5)
  #   updateSliderInput(session, "N", value = 2000)
  #   updateSliderInput(session, "input_noise_sd", value = .8)
  # })
  
  # Main center panel ----
  output$distPlot <- renderPlot({
    
    N <- input$N
    df <- vals$sim.df
    
    stats.df <- vals$stats.df
    
    p_t1 <- ggplot(df, aes(x=true, y=obs_t1, color=group)) + geom_point(size=3, alpha=.5) + theme(legend.position = "none") + 
      xlab("True score") + ylab("Observed score (Time 1)") +
      ggtitle(paste0("TRUE SCORE VS OBSERVED (1)\n",
                    "True Cohen's d: ", round(stats.df$d_true, 3), " (p=", round(stats.df$true_pval, 3), ")\n",
                     "Obs Cohen's d: ", round(stats.df$d_t1, 3), " (p=", round(stats.df$t1_pval, 3), ")\n",
                     "Corr between true and obs: ", round(stats.df$acc1, 3)))
    p_t2 <- ggplot(df, aes(x=true, y=obs_t2, color=group)) + geom_point(size=3, alpha=.5) + theme(legend.position = "none") + 
      xlab("True score") + ylab("Observed score (Time 2)") +
      ggtitle(paste0("TRUE SCORE VS OBSERVED (2)\n",
                    "True Cohen's d: ", round(stats.df$d_true, 3), " (p=", round(stats.df$true_pval, 3), ")\n",
                     "Obs Cohen's d: ", round(stats.df$d_t2, 3), " (p=", round(stats.df$t2_pval, 3), ")\n",
                     "Corr between true and obs: ", round(stats.df$acc2, 3)))
    p_rel <- ggplot(df, aes(x=obs_t1, y=obs_t2, color=group)) + geom_point(size=3, alpha=.5) + 
      xlab("Observed score (Time 1)") + ylab("Observed score (Time 2)") +
      ggtitle(paste0("TEST-RETEST RELIABILITY\n",
                    "Obs reliability for Gp 1: ", round(stats.df$c1, 3), "\n",
                     "Obs reliability for Gp 2: ", round(stats.df$c2, 3)))
    plot_grid(p_t1, p_t2, p_rel, ncol = 3, align = "hv")
    
  })
  
  output$distPlot100 <- renderPlot({
    
    df <- vals$stats100.df
    Nsim <- nrow(df)
    
    df_long <- data.frame(value=c(df$true_pval, df$t1_pval, df$t2_pval, df$d_true, df$d_t1, df$d_t2, df$c1, df$c2),
                          type=c(rep("pval", 3*Nsim), rep("d", 3*Nsim), rep("c", 2*Nsim)),
                          instance=c(rep(c(rep("Actual/True",Nsim), rep("Obs (Time 1)",Nsim), rep("Obs (Time 2)",Nsim)), 2), rep("Obs (Time 1)",Nsim), rep("Obs (Time 2)",Nsim))
                          )
     
    sigP_true <- sum(df$true_pval<.05)
    sigP_t1 <- sum(df$t1_pval<.05)
    sigP_t2 <- sum(df$t2_pval<.05)
    sigP_t1true <- sum(df$true_pval<.05 & df$t1_pval<.05)
    sigP_t2true <- sum(df$true_pval<.05 & df$t2_pval<.05)

    p_pval <- ggplot(df_long[which(df_long$type=="pval"), ], aes(x=instance, y=value)) + geom_violin() + ylab("p-value") + ylim(c(0,1)) +xlab("") +
      theme(axis.text.x = element_text(face="bold", color="black", size=12)) +
      ggtitle(paste0("GROUP DIFFERENCES DETECTED\n",
                     "True # p<.05: ", sigP_true, "\n",
                     "Obs # p<.05 (Time1): ", sigP_t1, " (", sigP_t1true, " true)\n",
                     "Obs # p<.05 (Time2): ", sigP_t2, " (", sigP_t2true, " true)\n"))
    p_d <- ggplot(df_long[which(df_long$type=="d"), ], aes(x=instance, y=value)) + geom_violin() + ylab("Cohen's d") + ylim(c(0,5)) +xlab("") +
      theme(axis.text.x = element_text(face="bold", color="black", size=12)) +
      ggtitle(paste0("EFFECT SIZES\n",
                     "Mean d (true): ", round(mean(df$d_true),3), "\n",
                     "Mean d (obs time 1): ", round(mean(df$d_t1),3), "\n",
                     "Mean d (obs time 2): ", round(mean(df$d_t2),3), "\n"))
    p_c <- ggplot(df_long[which(df_long$type=="c"), ], aes(x=instance, y=value)) + geom_violin() + ylab("Reliability") + ylim(c(0,1)) +xlab("") +
      theme(axis.text.x = element_text(face="bold", color="black", size=12)) +
      ggtitle(paste0("TEST-RETEST RELIABILITY\n",
                     "Mean reliability (group 1): ", round(mean(df$c1),3), "\n",
                     "Mean reliability (group 2): ", round(mean(df$c2),3), "\n"))
      scale_x_discrete(labels=c("Group 1", "Group 2"))

    plot_grid(p_pval, p_d, p_c, ncol = 3, align = "hv")
    
  })
  
  
  output$text_params <- renderText({
    paste("<br>",
          "<font size=2><p><b>Expected Group Difference:</b> Expected difference between groups, or size of effect in Group 1 where Group 2 has an effect of 0 (same for both time points).  This sets the means in the normal distributions the data are drawn from.</p>",
          "<p><b>True Population Variation:</b> Population variation on the measure of interest, assumed to be the same for both groups.  This sets the standard deviations in the normal distributions the data are drawn from.</p>",
          "<p><b>Sample Size:</b> Number of participants in each group</p>",
          "<p><b>Measurement Error:</b> Random amounts that observed scores differ from actual scores for both groups.</p></font>",
          sep="")
  })
  
  output$text <- renderText({
    
    df <- vals$stats.df

    paste("<h4>Results</h4>",
          "<p>Gp dif (Cohen's d) at time 1: <b>", round(df$d_t1, 3), "</b>.  p=", round(df$t1_pval, 3), "</b><br>",
          "Gp dif (Cohen's d) at time 2: <b>", round(df$d_t2, 3), "</b>.  p=", round(df$t2_pval, 3), "</b></p>",
          "(The true effect size is: <b>", round(df$d_true, 3), ")</p>",
          "<p>Reliability for Gp 1: <b>", round(df$c1, 3), "</b></p>",
          "<p>Reliability for Gp 2: <b>", round(df$c2, 3), "</b></p>",
          sep="")

  })
  
}


shinyApp(ui = ui, server = server)


