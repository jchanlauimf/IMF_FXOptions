# Structured Curriculum FX Options

rm(list=ls()) # Clean up memory
library(ggplot2) # Graphic library
library(lubridate) # Date manipulation library
library(dplyr) # Data manipulation library
library(RND)
source("auxFunctions.R")

my_wdir="D:/GitHub Projects/IMF_FXOptions/presentation"
setwd(my_wdir)

filename = "2018_IET_Options_data.csv"
data = read.csv(filename, header=TRUE)
data$Dates = mdy_hm(as.character(data$Dates))

# Plotting the spot exchange rate

p1 = ggplot(data, aes(Dates, Spot)) + geom_line(colour="blue")
p1 = p1 + geom_vline(xintercept=as.POSIXct(as.Date(("2016-06-22 UTC"))),
                     linetype="longdash")
p1

p2 = ggplot(data, aes(Dates, FWD3M)) + geom_line(colour="red")
p2 = p2 + geom_vline(xintercept=as.POSIXct(as.Date(("2016-06-22 UTC"))),
                     linetype="longdash")
p2


data$spreadRd = (data$ImpliedRd - data$Rd)*100

p3 = ggplot(data, aes(Dates, spreadRd)) + geom_line(colour="purple") +
        geom_hline(yintercept=0)
p3 = p3 + geom_vline(xintercept=as.POSIXct(as.Date(("2016-06-22 UTC"))),
                     linetype="longdash")
