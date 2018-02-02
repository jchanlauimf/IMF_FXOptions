# -------------------------------------------------------------------
# auxFunctions.R
#
# This file contains auxiliary functions accompanying the textbook:
#
# "Derivatives without Tears - A Guide for Non-Quant Professionals"
#
# and related set of lecture notes and courses. 
#
# Please do not distribute without authorization.
#
# Author: Jorge A. Chan-Lau
#
# First Version: February 1, 2017
# Last revision: December 4, 2017
# -------------------------------------------------------------------

library(grid)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

GKoption.premium = function(K,sigma,S,Tenor,fwd,rf,option_type)
{
  if (option_type =="c") {w=1}
  if (option_type =="p") {w=-1}
  d1 = log(fwd/K)+0.5*sigma*sigma*Tenor
  d1 = d1/(sigma*sqrt(Tenor))
  d2 = d1 - sigma*sqrt(Tenor)
  rd = log(fwd/S)/Tenor + rf
  premium = exp(-rd*Tenor)*(w*fwd*pnorm(w*d1) - w*K*pnorm(w*d2))
  return(premium)
}

