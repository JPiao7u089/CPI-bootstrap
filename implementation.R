rm(list=ls())
library(tidyverse)
library(survival)
library(dplyr)


##################################################################################################
################################      Example Codes    ###########################################
##################################################################################################

# load training data 
training_data = read.csv("training_data.csv")

# load testing data
testing_data = read.csv("testing_data.csv")

# load functions
source("myfunction.R")


#################################################################################################
## Conformal predictive confidence intervals
## working model: log-normal
## side: one
## alpha: 0.05
## output: coverage 
################################################################################################# 

conf.aft.fun(x=training_data[,c("x1","x2","x3","x4")],y=training_data$y,delta=training_data$delta,dat=training_data,
                              test.t=testing_data$t, test.x=testing_data[,c("x1","x2","x3","x4")],dat_test=testing_data,
                              alpha=0.05, side="one")

# [1] 0.9408186

#################################################################################################
## Conformal predictive confidence intervals
## working model: log-normal
## side: two
## alpha: 0.10
################################################################################################# 

conf.aft.fun(x=training_data[,c("x1","x2","x3","x4")],y=training_data$y,delta=training_data$delta,dat=training_data,
                            test.t=testing_data$t, test.x=testing_data[,c("x1","x2","x3","x4")],dat_test=testing_data,
                            alpha=0.10, side="two")
#[1] 0.903208

#################################################################################################
## Conformal predictive confidence intervals
## working model: Weibull
## side: two
## alpha: 0.10
################################################################################################# 

conf.weibull.fun(x=training_data[,c("x1","x2","x3","x4")],y=training_data$y,delta=training_data$delta,dat=training_data,
                 test.t=testing_data$t, test.x=testing_data[,c("x1","x2","x3","x4")],dat_test=testing_data,
                 alpha=0.10, side="two")
#[1] 0.903208

#################################################################################################
## conformal method
## working model: Cox
## side: two
## alpha: 0.10
################################################################################################# 
conf.cox.fun(x=training_data[,c("x1","x2","x3","x4")],y=training_data$y,delta=training_data$delta,dat=training_data,
             test.t=testing_data$t, test.x=testing_data[,c("x1","x2","x3","x4")],
              alpha=0.10, side="two")
#[1] 0.9021018