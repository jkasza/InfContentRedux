#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(reshape2)
library(plyr)
library(swCRTdesign)
library(matrixcalc)
library(scales)

source("InfContentInvestigations_functions.R", local=TRUE)


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
  titlePanel(h1("InfContentRedux",h3("Information content of SW designs with treatment effect heterogeneity and implementation periods"))),

   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
     
          sliderInput(
            "K", "Number of treatment sequences:",
            min = 3,
            max=20,
            step = 1,
            value = 3),
          
          radioButtons("impYN", label = ("Allow for implementation periods?"), 
                       choices = list("Yes" = 1, "No" = 2), selected = 2),
          conditionalPanel(
            condition = "input.impYN == '1'",
            uiOutput("sliderImp")#,
           #  sliderInput(
           #   "nIMP", "Number of periods required for implementation:",
           #   min = 1,
           #   max=10,
          #    step = 1,
          #    value = 1) 
             ),
        #sliderInput("T",
        #            "Number of periods, T:",
        #            min = 2,
        #            max=20,
        #            step = 1,
        #            value = 10),
        
        numericInput("m",
                     "Number of subjects in each cluster-period, m:",
                     min = 1,
                     max=1000,
                     step = 1,
                     value = 100),
        numericInput("icc_cc",
                     label= HTML(paste("&rho;",tags$sub("CC"),": ICC of subjects in the same control period", sep="")),
                     min = 0,
                     max=1,
                     step = 0.005,
                     value = 0.05),
        numericInput("icc_tt",
                     label= HTML(paste("&rho;",tags$sub("TT"),": ICC of subjects in the same intervention period", sep="")),
                     min = 0,
                     max = 1,
                     step = 0.005,
                     value = 0.05),
        numericInput("icc_ct",
                     label= HTML(paste("&rho;",tags$sub("CT"),": ICC of subjects, one in control, the other in intervention", sep="")),
                     min = 0,
                     max=1,
                     step = 0.005,
                     value = 0.05),
        numericInput("sige2",
                     "Residual variance",
                     min = 0,
                     step = 0.005,
                     value = 0.95),
        radioButtons("decayYN", label = ("Allow for decay in between-period ICC over time?"), 
                    choices = list("Yes" = 1, "No" = 2), selected = 2),
        conditionalPanel(
          condition = "input.decayYN == '1'",
          sliderInput(
            "decay", "Percentage decay per period:",
            min = 0,
            max=1,
            step = 0.001,
            value = 0.05) )
        ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel("Information content of cells",  
                   plotOutput("varplot"),
                   textOutput("ICcelltext")
                  ),
          tabPanel("Information content of sequences",  
                   plotOutput("varSEQplot"),
                   textOutput("ICseqtext")
                  ),
          tabPanel("Information content of periods",  
                   plotOutput("varPERplot"),
                   textOutput("ICpertext")
                  ),
          tabPanel("Trial schematic",  
                   plotOutput("desmatplot")
                  ),
          tabPanel("Reference and contact details",  
                   textOutput("refsdeets")
          )
          )
        )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$sliderImp <- renderUI({
    sliderInput(inputId = "nIMP", label = "Number of periods required for implementation:",
                min = 1, max = input$K-2, value = 1,  step = 1)
  })
  
  output$ICcelltext <- renderText({
    "The information content of sequence-period cells is displayed, where the information content of a cell
    is the variance of the treatment effect estimator if that cell is excluded divided by the variance if that 
    cell is included. There may be some cells which are entirely black:
    this implies that the model as specified cannot be fitted if this cell is excluded."
  })
  
  output$ICseqtext <- renderText({
    "The information content of treatment sequences is displayed, where the information content of a sequence
    is the variance of the treatment effect estimator if that sequence is excluded divided by the variance if that 
    sequence is included. There may be some sequences which are entirely black:
    this implies that the model as specified cannot be fitted if this sequence is excluded."
  })
  
  output$ICpertext <- renderText({
    "The information content of periods is displayed, where the information content of a period
    is the variance of the treatment effect estimator if that period is excluded divided by the variance if that 
    period is included. Gray period are those missing by design. There may be some periods which are entirely black:
    this implies that the model as specified cannot be fitted if this period is excluded."
  })
  
  output$refsdeets <- renderText({
    "This app accompanies the following paper:
    
    Kasza J, Taljaard M, Forbes AB. Information content of stepped-wedge designs when treatment
    effect heterogeneity and/or implementation periods are present. Submitted 2018.
    
    Please contact Jessica Kasza, jessica.kasza@monash.edu, with comments or questions.
    "
  })
  
  
  
  output$varplot <- renderPlot({
    
    #if(input$design == 1) Xdes <- SWdesmat(input$T)
    #if(input$design == 3) Xdes <- pllelbasedesmat2(input$T, input$K)
    #if(input$design == 4) Xdes <- crxodesmat2(input$T, input$K)
    #if(input$design == 5) Xdes <- SWdesmat2sandwich(input$T)
    Xdes <- SWdesmat(input$K + 1)
    if(input$impYN==1) {
      Xdes <- SWdesmatIMP(input$K + 1, input$nIMP)
    }
    #If there are any periods with only 2 treatment sequences, with 
    #differential exposure, cannot calculate the information content 
    #of either of these two sequences.
    #Need to flag that the cells in these periods MUST be included and
    #thus do not have an information content.
    Xdes_colsum <- colSums(Xdes, na.rm=TRUE)
    Xdes_nasums <- colSums(is.na(Xdes) == FALSE)
    
    
    #Generate the required variances from the ICCs
    sigu2 <- input$icc_cc*input$sige2/(1-input$icc_cc)
    combo <- (input$icc_tt*(sigu2 + input$sige2) - sigu2)/(1-input$icc_tt)
    siguv <- input$icc_ct*sqrt(sigu2 + combo + input$sige2)*sqrt(sigu2 + input$sige2) - sigu2
    sigv2 <- combo - 2*siguv
    
    myvar <- matrix(c(sigu2, siguv, siguv, sigv2), nrow=2)
    
    validate(
      need(is.positive.semi.definite(myvar), "Selected intra-cluster correlation values do not lead to a positive semi-definite variance matrix. Try a smaller value for rho_CT.")
    )
    
    
    varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
    
    if(input$decayYN == 2)  varmatall <- CRTVar_TEHet(Xdes, input$m, sigu2, sigv2, siguv, input$sige2)
    else if(input$decayYN == 1){
      myr = 1- input$decay
      varmatall <- CRTVar_TEHet_decay(Xdes, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)
    }
    #Code in for loops, but re-do this later
    for(i in 1:nrow(Xdes)) {
      for(j in 1:ncol(Xdes)) {
        if(is.na(Xdes[i,j])==TRUE)  varmat_excl[i,j] <- NA
        else if(is.na(Xdes[i,j])==FALSE) {
          if(Xdes_colsum[j] == 1 & Xdes_nasums[j]==2) varmat_excl[i,j] <- 2.2772
          else   {
            Xdesij <- Xdes
            Xdesij[i,j] <- NA
            if(input$decayYN == 2)  varmat_excl[i,j] <- CRTVar_TEHet(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2)/varmatall 
            else if(input$decayYN == 1)  varmat_excl[i,j] <- CRTVar_TEHet_decay(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)/varmatall
          }
        }
      }
    }
    
    varmat_excl<-round(varmat_excl, 4)
    melted_varmatexcl <- melt(varmat_excl)
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var1"] <- "Sequence"
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var2"] <- "Period"

    color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))
    if(sum(melted_varmatexcl$value==2.2772, na.rm=TRUE) > 0)    color_palette[length(table(melted_varmatexcl$value))]<- "#000000"
    
    T <- ncol(Xdes)
    K <- nrow(Xdes)
    
    ggplot(data = melted_varmatexcl, aes(x=Period, y=Sequence, fill = factor(value))) + 
      geom_tile( colour = "grey50") +
      scale_y_reverse(breaks=c(1:K)) +
      scale_x_continuous(breaks=c(1:T)) +
      theme(panel.grid.minor = element_blank()) +      
      geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 4) +
      scale_fill_manual(values = color_palette, breaks=levels(melted_varmatexcl$value)[seq(90, 150, by=5)]) 
    
  })
  
  
  output$varSEQplot <- renderPlot({
    #if(input$design == 1) Xdes <- SWdesmat(input$T)
    #if(input$design == 3) Xdes <- pllelbasedesmat2(input$T, input$K)
    #if(input$design == 4) Xdes <- crxodesmat2(input$T, input$K)
    #if(input$design == 5) Xdes <- SWdesmat2sandwich(input$T)
    Xdes <- SWdesmat(input$K + 1)
    if(input$impYN==1) {
      Xdes <- SWdesmatIMP(input$K + 1, input$nIMP)
    }
    
    #Generate the required variances from the ICCs
    sigu2 <- input$icc_cc*input$sige2/(1-input$icc_cc)
    combo <- (input$icc_tt*(sigu2 + input$sige2) - sigu2)/(1-input$icc_tt)
    siguv <- input$icc_ct*sqrt(sigu2 + combo + input$sige2)*sqrt(sigu2 + input$sige2) - sigu2
    sigv2 <- combo - 2*siguv
    
    myvar <- matrix(c(sigu2, siguv, siguv, sigv2), nrow=2)
    
    validate(
      need(is.positive.semi.definite(myvar), "Selected intra-cluster correlation values do not lead to a positive semi-definite variance matrix. Try a smaller value for rho_CT.")
    )
    
    
    varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
    
    if(input$decayYN == 2)  varmatall <- CRTVar_TEHet(Xdes, input$m, sigu2, sigv2, siguv, input$sige2)
    else if(input$decayYN == 1){
      myr = 1- input$decay
      varmatall <- CRTVar_TEHet_decay(Xdes, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)
    }
    
    #Code in for loops, but re-do this later
    if(input$impYN==2) {
    
       for(i in 1:nrow(Xdes)) {
          Xdesij <- Xdes
          Xdesij[i,] <- NA
          if(input$decayYN == 2)  varmat_excl[i,] <- CRTVar_TEHet(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2)/varmatall 
          else if(input$decayYN == 1)  varmat_excl[i,] <- CRTVar_TEHet_decay(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)/varmatall
          }
    }
    else if(input$impYN==1) {
        if(input$nIMP == (input$K - 2)) {
          #Then cannot remove the outer sequences
          varmat_excl[1,] <- 34
          varmat_excl[nrow(Xdes),] <- 34
          for(i in 2:(nrow(Xdes)-1)) {
            Xdesij <- Xdes
            Xdesij[i,] <- NA
            if(input$decayYN == 2)  varmat_excl[i,] <- CRTVar_TEHet(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2)/varmatall 
            else if(input$decayYN == 1)  varmat_excl[i,] <- CRTVar_TEHet_decay(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)/varmatall          }
        }
        else if(input$nIMP < (input$K - 2)){
          
          for(i in 1:nrow(Xdes)) {
            Xdesij <- Xdes
            Xdesij[i,] <- NA
            if(input$decayYN == 2)  varmat_excl[i,] <- CRTVar_TEHet(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2)/varmatall 
            else if(input$decayYN == 1)  varmat_excl[i,] <- CRTVar_TEHet_decay(Xdesij, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)/varmatall          }
          
        }
    }
    
    varmat_excl<-round(varmat_excl, 4)
    melted_varmatexcl <- melt(varmat_excl)
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var1"] <- "Sequence"
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var2"] <- "Period"

    color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))
    if(sum(melted_varmatexcl$value==34, na.rm=TRUE) > 0)    color_palette[length(table(melted_varmatexcl$value))]<- "#000000"
    
    T <- ncol(Xdes)
    K <- nrow(Xdes)
    
    ggplot(data = melted_varmatexcl, aes(x=Period, y=Sequence, fill = factor(value))) + 
      geom_tile( colour = "grey50") +
      scale_y_reverse(breaks=c(1:K)) +
      scale_x_continuous(breaks=c(1:T)) +
      theme(panel.grid.minor = element_blank()) +         
      geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 4) +
      scale_fill_manual(values = color_palette, breaks=levels(melted_varmatexcl$value)[seq(90, 150, by=5)])
  })
  
  
  output$varPERplot <- renderPlot({
    #if(input$design == 1) Xdes <- SWdesmat(input$T)
    #if(input$design == 3) Xdes <- pllelbasedesmat2(input$T, input$K)
    #if(input$design == 4) Xdes <- crxodesmat2(input$T, input$K)
    #if(input$design == 5) Xdes <- SWdesmat2sandwich(input$T)
    Xdes <- SWdesmat(input$K + 1)
    if(input$impYN==1) {
      Xdes <- SWdesmatIMP(input$K + 1, input$nIMP)
    }
    
    #Generate the required variances from the ICCs
    sigu2 <- input$icc_cc*input$sige2/(1-input$icc_cc)
    combo <- (input$icc_tt*(sigu2 + input$sige2) - sigu2)/(1-input$icc_tt)
    siguv <- input$icc_ct*sqrt(sigu2 + combo + input$sige2)*sqrt(sigu2 + input$sige2) - sigu2
    sigv2 <- combo - 2*siguv
    
    myvar <- matrix(c(sigu2, siguv, siguv, sigv2), nrow=2)
    
    validate(
      need(is.positive.semi.definite(myvar), "Selected intra-cluster correlation values do not lead to a positive semi-definite variance matrix. Try a smaller value for rho_CT.")
    )
    
    
    varmat_excl<-matrix(data=NA, nrow=nrow(Xdes), ncol=ncol(Xdes))
    
    if(input$decayYN == 2)  varmatall <- CRTVar_TEHet(Xdes, input$m, sigu2, sigv2, siguv, input$sige2)
    else if(input$decayYN == 1){
      myr = 1- input$decay
      varmatall <- CRTVar_TEHet_decay(Xdes, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)
    }    
    #If there is only one period during which clusters are exposed to 
    #each condition, then this period cannot be omitted from the study
    Xdes0 <- (colSums(Xdes==0, na.rm=TRUE) ==1)
    Xdes1 <- (colSums(Xdes==1, na.rm=TRUE) ==1)
    Xdes01 <- colSums(rbind(Xdes0, Xdes1))
    
        
    #Code in for loops, but re-do this later
    if(input$decayYN == 2) {
       for(j in 1:ncol(Xdes)) {      
         if(Xdes01[j] != 2)    varmat_excl[,j] <- CRTVar_TEHet_ocol(Xdes, j, input$m, sigu2, sigv2, siguv, input$sige2)/varmatall 
         else if (Xdes01[j] == 2) varmat_excl[,j] <- 34
      }
    }
    else if(input$decayYN == 1){
      for(j in 1:ncol(Xdes)) {      
        if(Xdes01[j] != 2)    varmat_excl[,j] <- CRTVar_TEHet_ocol_decay(Xdes, j, input$m, sigu2, sigv2, siguv, input$sige2, type=1, r=myr)/varmatall 
        else if (Xdes01[j] == 2) varmat_excl[,j] <- 34
      }
      
    }    
    varmat_excl<-round(varmat_excl, 4)
    melted_varmatexcl <- melt(varmat_excl)
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var1"] <- "Sequence"
    names(melted_varmatexcl)[names(melted_varmatexcl)=="Var2"] <- "Period"

    color_palette <-colorRampPalette(c( "yellow", "red"))(length(table(varmat_excl)))
    if(sum(melted_varmatexcl$value==34, na.rm=TRUE) > 0)    color_palette[length(table(melted_varmatexcl$value))]<- "#000000"
    
    T <- ncol(Xdes)
    K <- nrow(Xdes)
    
    ggplot(data = melted_varmatexcl, aes(x=Period, y=Sequence, fill = factor(value))) + 
      geom_tile( colour = "grey50") +
      scale_y_reverse(breaks=c(1:K)) +
      scale_x_continuous(breaks=c(1:T)) +
      theme(panel.grid.minor = element_blank()) +           
      geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 4) +
      scale_fill_manual(values = color_palette, breaks=levels(melted_varmatexcl$value)[seq(90, 150, by=5)])
  })
   
  
  output$desmatplot <- renderPlot({
    Xdes <- SWdesmat(input$K + 1)
    if(input$impYN==1) {
      Xdes <- SWdesmatIMP(input$K + 1, input$nIMP)
    }
    
    melted_myx <- melt(Xdes)
    names(melted_myx)[names(melted_myx)=="Var1"] <- "Sequence"
    names(melted_myx)[names(melted_myx)=="Var2"] <- "Period"
    
    color_palette <-colorRampPalette(c( "white", "grey"))(2)
    
    K <- nrow(Xdes)
    T <- ncol(Xdes)
    
    ggplot(data =melted_myx, aes(x=Period, y=Sequence, fill = factor(value))) + 
      geom_tile( colour = "grey50") +
      scale_y_reverse(breaks=c(1:K)) +
      scale_x_continuous(breaks=c(1:T)) +
      theme(panel.grid.minor = element_blank()) +           
      geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 10) +
      scale_fill_manual(values = color_palette) +  theme(legend.position="none")
    
  })  
}

# Run the application 
shinyApp(ui = ui, server = server)

