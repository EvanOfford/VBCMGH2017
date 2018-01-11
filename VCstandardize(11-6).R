

require(nplr)
library(nplr) #this package is for running N-parameter logistic regression analysis

directory <- "C:/Users/eso7/Documents/R/VCDATA/Outlier_Tests"#Direct a path to where the .csv files are stored. The last file is your folder, ex. 'Outlier Tests'
analyze <- function(filename = NULL) {
dat <- read.csv(file = filename, header = FALSE, fileEncoding="UTF-8-BOM")
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
custom_serumDfac <- c()#Input custom vibriocidal dilutions and change the variable on line 27
mAb <- c(2500, 1250, 625, 312.5, 156.25, 78.125, 39.0625, 19.53125, 9.765625, 4.8828125, 2.44140625)
custom_mAb <- c(1000,500,250,125,62.5,31.25,15.625,7.8125,3.90625,1.953125,0.9765625)#Input custom monoclonal antibody dilutions and change the variable on line 28
serumlog <- log10(serumDfac)
mAblog <- log10(custom_mAb)

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

PlateSTD <- as.numeric(c(STD_VB_TITER_P1S1, STD_VB_TITER_P1S2, STD_VB_TITER_P1S3))
return(PlateSTD)

}

filename <- list.files(path = directory, pattern = ".csv", full.names = TRUE, include.dirs = F)
filenames <- list.files(path = directory, pattern = ".csv", full.names = F, include.dirs = F)
filenames <- gsub(".csv", " ", filenames)

Test <- matrix(data = (lapply(filename, analyze)), nrow = length(filenames), ncol = 3)
Test[1:length(filenames),2] <- ""
Test[1:length(filenames),3] <- ""

rownames(Test) <- paste0(filenames)
colnames(Test) <- as.character(c("Sample 1", "Sample 2", "Sample 3"))

setwd("C:/Users/eso7/Documents/R/Results_Folder")#Direct a path to the folder where you would like the results file to be stored. The last file is your folder. ex. 'Results_Folder'.
write.csv(Test, file = paste0(Sys.Date()," ", "Standardized_Data.csv"))#Name your file (be sure to end with .csv).




