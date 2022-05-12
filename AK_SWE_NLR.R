#This model uses a nonlinear regression to predict SWE
#The code used to generate and fit the model was written by Eric Sproles, the rest was written by Sam Neitlich

#Set working directory
setwd("~/Documents/GPHY484_491/491Proj") #Change to your wd


#Install and library packages
install.packages('tidyverse')
install.packages('hydroGOF')
install.packages('zeallot')

library(tidyverse)
library(caret)
library('hydroGOF')
library(ggplot2)
library(stringr)
library(raster)
library(zeallot)


## read in files from snotel sites to develop parameters with a subset of the sites
BASE1 <-read.csv(file.choose()) #Read in nlr_inputs_base.csv here

#Define functions based on the equations in Leibowitz et al.
S_EXP  <- function(PRISM_T,S1,S2) (exp(-(PRISM_T + S1)/S2))
C_SNOW <- function(PRISM_T,S1,S2) (1 - (1/(1 + S_EXP(PRISM_T,S1,S2) )))
M_EXP  <- function(PRISM_T,M1,M2) (exp(-(PRISM_T + M1)/M2))
C_MELT <- function(PRISM_T,M1,M2) (1/(1 + M_EXP(PRISM_T,M1,M2)))

#Fit a reverse sigmoid nonlinear regression to the data
baseswe.mod <- nls(SNOTEL_SWE ~ SURPLUS * C_SNOW(PRISM_T,S1,S2) 
                   + PREV_SWE * (1 - C_MELT(PRISM_T,M1,M2)),
                   data=BASE1,start=list(S1 = -4, S2 = 2, M1 = -3, M2 = 1),trace = TRUE)

#S1 = -4, S2 = 5, M1 = -10, M2 = 0.1 are the original parameters used in the Bristol Bay model
#I updated them to S1 = -4, S2 = 2, M1 = -3, M2 = 1 for this model, but they are locationally dependant

#Summarize
summary(baseswe.mod)
#This will show us all of the parmater predictions that were optimized in this model

#Predict SWE to plot--note that we're predicting on the same dataset we trained the model on here 
BASE1['Predicted_SWE'] = predict(baseswe.mod)

#Generate linear model of predicted vs. observed SWE and summarize model
attach(BASE1)
lmod = lm(Predicted_SWE ~ SNOTEL_SWE)
summary(lmod)



#Plot the linear model in greater detail and fit regression line
ggplot(BASE1, aes(SNOTEL_SWE, Predicted_SWE))+
  geom_point()+
  labs(title="Predicted vs. Observed SWE", y= "Predicted SWE (mm)", x = "Observed SWE (mm)")+
  geom_smooth(method="lm")+
  theme_classic()+
  annotate("text", x=c(73,27,50), y=c(215,200,180), 
           label=c("Pred. SWE = 0.883 (Obs. SWE) + 1.731", 
                    expression(R^"2": 0.9698),
                   "Adjusted p-value: 2.2e-16"))




## This is the same as above but uses all of the SNOTEL sites in the model

COMB1 <-read.csv(file.choose()) #Read in file containing all sites (nlr_inputs_comb.csv)

#Define functions
S_EXP  <- function(PRISM_T,S1,S2) (exp(-(PRISM_T + S1)/S2))
C_SNOW <- function(PRISM_T,S1,S2) (1 - (1/(1 + S_EXP(PRISM_T,S1,S2) )))
M_EXP  <- function(PRISM_T,M1,M2) (exp(-(PRISM_T + M1)/M2))
C_MELT <- function(PRISM_T,M1,M2) (1/(1 + M_EXP(PRISM_T,M1,M2)))

#Fit regression
combswe.mod <- nls(SNOTEL_SWE ~ SURPLUS * C_SNOW(PRISM_T,S1,S2) 
                   + PREV_SWE * (1 - C_MELT(PRISM_T,M1,M2)),
                   data=COMB1,start=list(S1 = -4, S2 = 2, M1 = -3, M2 = 1),trace = TRUE)


#Summarize model and use to predict outputs
summary(combswe.mod)
COMB1['Predicted_SWE'] = predict(combswe.mod)

#Fit linear model and display with ggplot2-- same as above
attach(COMB1)
lmod = lm(Predicted ~ SNOTEL_SWE)
summary(lmod)

ggplot(COMB1, aes(SNOTEL_SWE, Predicted_SWE))+
  geom_point()+
  labs(title="Predicted vs. Observed SWE", y= "Predicted SWE (mm)", x = "Observed SWE (mm)")+
  geom_smooth(method="lm")+
  theme_classic()+
  annotate("text", x=c(73,27,50), y=c(215,200,180), 
           label=c("Pred. SWE = 0.883 (Obs. SWE) + 1.731", 
                   expression(R^"2": 0.9698),
                   "Adjusted p-value: 2.2e-16"))




#Create a random sample of the data and split into separate training and testing datasets
#This will provide a more robust analysis of model performance than what we did above. 
random_sample <- createDataPartition(COMB1$SNOTEL_SWE,
                                     p = 0.8, list = FALSE)
#Define training and testing datasets
training_dataset  <- COMB1[random_sample, ]
testing_dataset <- COMB1[-random_sample, ]

#Fit the model, same as above
S_EXP  <- function(PRISM_T,S1,S2) (exp(-(PRISM_T + S1)/S2))
C_SNOW <- function(PRISM_T,S1,S2) (1 - (1/(1 + S_EXP(PRISM_T,S1,S2) )))
M_EXP  <- function(PRISM_T,M1,M2) (exp(-(PRISM_T + M1)/M2))
C_MELT <- function(PRISM_T,M1,M2) (1/(1 + M_EXP(PRISM_T,M1,M2)))

#Use training dataset as input this time
training.mod <- nls(SNOTEL_SWE ~ SURPLUS * C_SNOW(PRISM_T,S1,S2) 
                   + PREV_SWE * (1 - C_MELT(PRISM_T,M1,M2)),
                   data=training_dataset,start=list(S1 = -4, S2 = 2, M1 = -3, M2 = 1),trace = TRUE)

#Predict values and calculate four statistics
predictions <- predict(training.mod, testing_dataset)

data.frame( R2 = R2(predictions, testing_dataset $ SNOTEL_SWE),
            RMSE = RMSE(predictions, testing_dataset $ SNOTEL_SWE),
            MAE = MAE(predictions, testing_dataset $ SNOTEL_SWE),
            NSE = NSE(predictions, testing_dataset $ SNOTEL_SWE))

#These statistics indicate that the model is performing well

#-------------------------------------------------

coefficients <- summary(combswe.mod)$coefficients  # Extract coefficients in matrix
c(S1, S2, M1, M2) %<-% coefficients[c(1, 2, 3, 4)] #Set the coefficients equal to variables

#Import the zip file that we downloaded in python here
zipF<-file.choose() #File should be in downloads section

#Unzip the file to a directory (working directory is preferable)
outDir<-"/Users/samneitlich/Documents/GPHY484_491/491Proj/" # Define the folder where the zip file should be unzipped to 
unzip(zipF,exdir=outDir)




list <- list(10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9) #List months by water year date

for (i in list) { #Iterate through every month
  
  #Put in a leading 0 and determine the file name to import for that month
  s <- str_pad(i, 2, pad = "0")
  fname = paste("nlr_input", s, ".csv", sep = "")
  
  #Read the csv and save to dataframe
  input <- read.csv(file = fname, header = TRUE)
  #nrow(input)
  
  #Set the previous SWE for October to 0, otherwise, use the predicted SWE from the previous month
  #This is why we sorted the months in order of the water year
  if (i==10) {
    input['PREV_SWE'] <- 0
  } else {
    input['PREV_SWE'] <- toRaster['SWE']
  }

  
  #Create a new dataframe with only the values used in prediction, this is necessary for the predict function to work
  new <- data.frame(SURPLUS=c(input$SURPLUS),PRISM_T=c(input$PRISM_T), PREV_SWE=c(input$PREV_SWE))
  
  #Predict the SWE for each location in the watershed
  pred <- predict(baseswe.mod, newdata = new)
  
  #Save predictions to column of dataframe
  input['PredictedSWE'] <- pred
  
  #Change all of the predicted SWE values that are below 0 to 0
  input$PredictedSWE <- ifelse(input$PredictedSWE >0, input$PredictedSWE, 0)
  #Use the same logic to change any values above 300 to 300, as this is roughly that max SWE that we would expect per the SNOTEL site data
  input$PredictedSWE <- ifelse(input$PredictedSWE <300, input$PredictedSWE, 300)
  
  #Create a new dataframe with only lat, lon, and predicted SWE
  toRaster <- data.frame(X=c(input$X),Y=c(input$Y), SWE=c(input$PredictedSWE))
  
  rst <- rasterFromXYZ(toRaster, crs = 4269) #Create raster from dataframe, crs = NAD 83
  rst
  
  title = paste("Distribution of SWE values, month = ", i, sep="") #Set title
  
  
  #Plot a histogram of values in the raster
  hist(rst, main=title,
       xlab = 'Predicted SWE (mm)',
       ylab = 'Count',
       col= "blue",
       maxpixels=22000000)

  plot(rst, main=title) #Plot the raster

}

#Now we can see histograms and rasters of the SWE values in this watershed for every month
#These files can easily be exported or saved. 

#Some of the lower values are getting overpredicted, and I keep rounding them down to 300. 
#While this keeps the error from propegating too far, it's not a perfect solution. 
#That said, the model performs pretty well throughout most of the watershed. 

