library(rsconnect)
library(shiny)
library(nplr)#this package is for running N-parameter logistic regression analysis
library(shiny)
library(shinythemes)



ui <- navbarPage(theme = shinytheme("lumen"),
  "Vibriocidal Standardization with a Monoclonal Antibody",
                 
                 tabPanel("Application", div(img(src = "mgh1811.png", height = 200, width = 200), style="text-align: right;"),
                          
                          
                          mainPanel(div(style="text-align: left; margin-top: -13em", 
                                        
                                        fileInput('file1', 'Choose a .CSV file to Statistically Analyze.', multiple = F,
                                                  accept = c('text/csv',
                                                             'text/comma-separated-values',
                                                             'text/tab-separated-values',
                                                             'text/plain','.csv','.tsv')))),
                          
                          
                          uiOutput('monoclonalUI'), br(),
                          uiOutput('buttonsUI'), br(),
                          uiOutput('downloadUI'), br()
                        
                     
                 ),
                 
                 tabPanel("Instructions for First Time Use",
                          div(img(src = "VBC1.png", height = 1000, width = 880), style="text-align: left; margin-top: -5em"),
                          div(img(src = "VBC2.png", height = 1000, width = 880), style="text-align: left; margin-top: -13em"),
                          div(img(src = "Plate.png", height = 300, width = 800), style="text-align: left; margin-top: -18em; margin-bottom: 2em"))
)

server = shinyServer(function(input, output) {
  
  controlVar <- reactiveValues(fileReady = FALSE, monoclonalReady = FALSE, interpReady = FALSE)
  
  dat <- NULL
  
  observeEvent(input$file1, {
    controlVar$fileReady <- FALSE
    if (is.null(input$file1))
      return()
    inFile <- input$file1
    dat <<- read.csv(inFile$datapath, header = F, fileEncoding="UTF-8-BOM")
    if(!is.data.frame(dat))
      return()
    controlVar$fileReady <- TRUE
  })
  output$monoclonalUI <- renderUI({
    if (controlVar$fileReady)
      div(
        numericInput('monoclonal', "Input your starting monoclonal antibody concentration in ng/mL. Default is 2500 ng/mL.", 2500), style="text-align: left; margin-top: -6em")
    
  })
  
  output$buttonsUI <- renderUI({
    if (controlVar$monoclonalReady && controlVar$fileReady)
      div(
        
        
        actionButton('Interpolate','Interpolate', style="text-align: left; margin-top: 0em"), style="color: #FF0000")
        
      
  })
  
  
  output$downloadUI <- renderUI({
    if (controlVar$interpReady && controlVar$monoclonalReady && controlVar$fileReady)
      div(
        downloadButton('downloadData', "Download the Data!"), style="text-align: left; margin-top: 0em")
  })
  
  
  
  observeEvent(input$monoclonal, {
    controlVar$monoclonalReady <- FALSE
    
    a <- (as.integer(input$monoclonal))/2
    b <- a/2
    c <- b/2
    d <- c/2
    e <- d/2
    f <- e/2
    g <- f/2
    h <- g/2
    i <- h/2
    j <- i/2
    
    mAb <<- c(S,a,b,c,d,e,f,g,h,i,j)
    
    controlVar$monoclonalReady <- TRUE
    
  })
  
  
  
  observeEvent(input$Interpolate, {
    controlVar$interpReady <- FALSE
    
    p1 <- array(data = dat)
    GC <- as.numeric(colMeans(p1[1:4,1:1]))
    negAve <- as.numeric(colMeans(p1[5:8,1:1]))
    fiftyGC1 <- (GC - negAve)/2
    mAbA2p1 <- as.numeric(colMeans(p1[7:8,2:12]))
    s1p1 <- as.numeric(colMeans(p1[1:2,2:12]))
    s2p1 <- as.numeric(colMeans(p1[3:4,2:12]))
    s3p1 <- as.numeric(colMeans(p1[5:6,2:12]))
    negVec <- rep(negAve, 11)
    mAbA2p1 <- mAbA2p1 - negVec
    S1p1 <- s1p1 - negVec
    S2p1 <- s2p1 - negVec
    S3p1 <- s3p1 - negVec
    
    serumDfac <- c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240)
    serumlog <- log10(serumDfac)
    mAblog <- log10(mAb)
    
    nplrP1S1 <- nplr(serumlog, S1p1, useLog = FALSE, LPweight = 0, npars = 4,
                     method = c("res"), silent = TRUE)
    nplrP1S2 <- nplr(serumlog, S2p1, useLog = FALSE, LPweight = 0, npars = 4,
                     method = c("res"), silent = TRUE)
    nplrP1S3 <- nplr(serumlog, S3p1, useLog = FALSE, LPweight = 0, npars = 4,
                     method = c("res"), silent = TRUE)
    
    nplrmAbP1 <- nplr(mAblog, mAbA2p1, useLog = FALSE, LPweight = 0, npars = 4,
                      method = c("res"), silent = TRUE)
    
    ADJ_VB_Titer_P1S1 <- getEstimates(nplrP1S1, target = fiftyGC1, B = 1e4, conf.level = .95)
    ADJ_VB_Titer_P1S2 <- getEstimates(nplrP1S2, target = fiftyGC1, B = 1e4, conf.level = .95)
    ADJ_VB_Titer_P1S3 <- getEstimates(nplrP1S3, target = fiftyGC1, B = 1e4, conf.level = .95)
    
    
    ADJ_mAb_P1 <- getEstimates(nplrmAbP1, target = fiftyGC1, B = 1e4, conf.level = .95)
    
    if(ADJ_VB_Titer_P1S1$x < 0.698970004) {
      ADJ_VB_Titer_P1S1$x <- 0.698970004 
    }
    if(ADJ_VB_Titer_P1S2$x < 0.698970004) {
      ADJ_VB_Titer_P1S2$x <- 0.698970004 
    }
    if(ADJ_VB_Titer_P1S3$x < 0.698970004) {
      ADJ_VB_Titer_P1S3$x <- 0.698970004 
    }
    
    if(ADJ_VB_Titer_P1S1$x == (1+1.6864291e-11) && s1p1[1] > .10) ADJ_VB_Titer_P1S1$x <- 0.69897000
    if(ADJ_VB_Titer_P1S2$x == (1+1.6864291e-11) && s2p1[1] > .10) ADJ_VB_Titer_P1S2$x <- 0.69897000
    if(ADJ_VB_Titer_P1S3$x == (1+1.6864291e-11) && s3p1[1] > .10) ADJ_VB_Titer_P1S3$x <- 0.69897000
    
    if(is.na(ADJ_VB_Titer_P1S1$x) == T && s1p1[1] > 0.10) {
      ADJ_VB_Titer_P1S1$x <- 0.69897000
    }
    if(is.na(ADJ_VB_Titer_P1S2$x) == T && s2p1[1] > 0.10) {
      ADJ_VB_Titer_P1S2$x <- 0.69897000
    }
    if(is.na(ADJ_VB_Titer_P1S3$x) == T && s3p1[1] > 0.10) {
      ADJ_VB_Titer_P1S3$x <- 0.69897000
    }
    
    if(ADJ_VB_Titer_P1S1$x > 2 && s1p1[1] > 0.10) {
      ADJ_VB_Titer_P1S1$x <- 0.69897000
    }
    if(ADJ_VB_Titer_P1S2$x > 2 && s2p1[1] > 0.10) {
      ADJ_VB_Titer_P1S2$x <- 0.69897000
    }
    if(ADJ_VB_Titer_P1S3$x > 2 && s3p1[1] > 0.10) {
      ADJ_VB_Titer_P1S3$x <- 0.69897000
    }
    
    if(is.na(ADJ_VB_Titer_P1S1$x) == T && s1p1[11] < 0.14) {
      ADJ_VB_Titer_P1S1$x <- 4.010299957
    }
    if(is.na(ADJ_VB_Titer_P1S2$x) == T && s2p1[11] < 0.14) {
      ADJ_VB_Titer_P1S2$x <- 4.010299957
    }
    if(is.na(ADJ_VB_Titer_P1S3$x) == T && s3p1[11] < 0.14) {
      ADJ_VB_Titer_P1S3$x <- 4.010299957
    }
    
    if(ADJ_VB_Titer_P1S1$x == (1+1.6864291e-11) && s1p1[11] < .14) ADJ_VB_Titer_P1S1$x <- 4.010299957
    if(ADJ_VB_Titer_P1S2$x == (1+1.6864291e-11) && s2p1[11] < .14) ADJ_VB_Titer_P1S2$x <- 4.010299957
    if(ADJ_VB_Titer_P1S3$x == (1+1.6864291e-11) && s3p1[11] < .14) ADJ_VB_Titer_P1S3$x <- 4.010299957
    
    STD_VB_TITER_P1S1 <- log2(10^ADJ_VB_Titer_P1S1$x) + log2(10^ADJ_mAb_P1$x)
    STD_VB_TITER_P1S2 <- log2(10^ADJ_VB_Titer_P1S2$x) + log2(10^ADJ_mAb_P1$x)
    STD_VB_TITER_P1S3 <- log2(10^ADJ_VB_Titer_P1S3$x) + log2(10^ADJ_mAb_P1$x)
    
    dat <<- matrix(data = (c(STD_VB_TITER_P1S1, STD_VB_TITER_P1S2, STD_VB_TITER_P1S3)), nrow = 1, ncol = 3, byrow=T)
    
    
    rownames(dat, prefix = paste0(input$file1))
    colnames(dat) <- c("Sample 1", "Sample 2", "Sample 3")
    
    controlVar$interpReady <- TRUE  
    
  })
  
  output$downloadData <- downloadHandler(
    filename = paste0("Standardized on ", Sys.Date(), ",", input$file1),
    content = function(file) {
      write.csv(dat, file, na="")
      
    }
  )
  
})


shinyApp(ui=ui, server=server)