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

library(fOptions)
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

greekOption = function(spot,option_par){
  # Function greekOption returns the desired greek given
  # the spot price of the underlying and the option characteristics
  
  r = as.numeric(option_par[1])
  b = as.numeric(option_par[2])
  K = as.numeric(option_par[3])
  H = as.numeric(option_par[4])
  t = as.numeric(option_par[5])
  sigma = as.numeric(option_par[6])
  option_type = option_par[7]
  greek_type = option_par[8]
  optionGreeks =  GBSCharacteristics(TypeFlag = option_type, S=spot, X = K, Time = H-t, 
                                     r =r, b = b, sigma = sigma)
  switch(greek_type,
         premium= optionGreeks$premium,
         delta  = optionGreeks$delta,
         theta  = optionGreeks$theta,
         vega   = optionGreeks$vega,
         rho    = optionGreeks$rho,
         lambda = optionGreeks$lambda,
         gamma  = optionGreeks$gamma
  )
}

BSMStrikeFromDelta = function(delta, sigma, option_par){
  
  # function calculates the strike corresponding to a
  # combination of delta and volatility
  
  S0 = as.numeric(option_par[1])
  H  = as.numeric(option_par[2])
  r  = as.numeric(option_par[3])
  b  = as.numeric(option_par[4])
  
  aux = qnorm(delta*exp(b*H))*sigma*sqrt(H)
  aux = aux -(r-b+1/2*sigma*sigma)*H
  K = S0/exp(aux)
  return(K);
}

functionRND = function(data.Delta, data.Ivol, data.other){
  
  # Read option characteristics
  
  S = data.other[1]
  H = data.other[2]
  r = data.other[3]
  b = data.other[4]
  t = data.other[5]
  
  # Fit second degree polynomial to the data
  fit.Vol = lm(data.Ivol~ poly(data.Delta,2,raw=TRUE))
  
  # Use fitted polynommial to interpolate Delta-Vol Curve
  delta.simple = seq(from=0.01, to = 0.99, by=0.005)
  delta.square = delta.simple*delta.simple
  delta.interc = rep(1,length(delta.simple))
  
  X = cbind(delta.interc, delta.simple, delta.square)
  iVolInterpol = t(t(X)*fit.Vol$coefficients)
  iVolInterpol = rowSums(iVolInterpol)
  
  # Check fit of the data visually

  
  # Recover the strike prices from the Deltas
  option_par =c(S,H,r,b)
  nstrikes = length(delta.simple)
  call.strike = rep(NA,nstrikes)
  for (i in 1:nstrikes){
    call.strike[i] = BSMStrikeFromDelta(delta.simple[i],
                                        iVolInterpol[i],
                                        option_par)
  }
  
  #plot(call.strike, iVolInterpol*100, xlab = "Strike", 
  #     ylab="Implied Volatility", type = "l")
  
  # Find call premium for each strike, required by RND
  option_type = "c"
  call.premium = rep(NA,nstrikes)
  greek_type = "premium"
  for (i in 1:nstrikes){
    option_par = c(r,b,call.strike[i],H,t,iVolInterpol[i],option_type,greek_type)
    call.premium[i] = greekOption(S,option_par)
  }
  
  # Find put premium for each strike, required by RND
  option_type = "p"
  put.premium = rep(NA,nstrikes)
  for (i in 1:nstrikes){
    option_par = c(r,b,call.strike[i],H,t,iVolInterpol[i],option_type,greek_type)
    put.premium[i] = greekOption(S,option_par)
  }
  
  call.strike  = sort(call.strike, decreasing=F)    # call strikes in ascending order
  call.premium = sort(call.premium, decreasing=T)   # call premium, in decreasing order
  put.strike   = call.strike                        # put strikes, same as call strikes
  put.premium  = sort(put.premium, decreasing=F)    # put premium in ascending order
  
  # Generalized Beta Distribution
  gb.rnd = extract.gb.density(r=r, te=H, y=b, s0=S, 
                              market.calls=call.premium, call.strikes = call.strike, call.weights =1,
                              market.puts = put.premium, put.strikes = put.strike, put.weights = 1,  
                              lambda=1, hessian.flag=F)
  
  # Mixed Normal Distribution
  mln.rnd = extract.mln.density(r=r, te=H, y=b, s0=S, 
                                market.calls=call.premium, call.strikes = call.strike, call.weights =1,
                                market.puts = put.premium, put.strikes = put.strike, put.weights = 1,  
                                lambda=1, hessian.flag=F, cl=list(maxit=10000))
  
  # EdgeWorth Expansion
  ew.rnd = extract.ew.density(initial.values = rep(NA,2), r=r, y = b, te=H, s0=S, 
                              market.calls=call.premium, call.strikes = call.strike, call.weights =1, 
                              lambda=1, hessian.flag=F, cl = list(maxit=10000))
  
  # Shimko Numerical Method
  shm.rnd = extract.shimko.density(market.calls=call.premium, call.strikes = call.strike, 
                                   r=r, y=b, t=H, s0=S, lower=-10, upper= 20)
  
  # Black-Scholes-Merton lognormal
  bsm.rnd = extract.bsm.density(initial.values = rep(NA,2), r=r, y = b, te=H, s0=S, 
                                market.calls=call.premium, call.strikes = call.strike, call.weights =1,
                                market.puts = put.premium, put.strikes = put.strike, put.weights = 1,  
                                lambda=1, hessian.flag=F, cl = list(maxit=10000))
  
  # Calculate the values of the risk-neutral distributions for the range of strikes
  K = seq(0.8*min(call.strike), 1.2*max(call.strike), 0.01)
  
  gb.points  = dgb(K,gb.rnd$a, gb.rnd$b, gb.rnd$v, gb.rnd$w)
  mln.points = dmln(K, mln.rnd$alpha.1, mln.rnd$meanlog.1, mln.rnd$meanlog.2,
                    mln.rnd$sdlog.1, mln.rnd$sdlog.2)
  ew.points  = dew(K,r,b,H,S,ew.rnd$sigma, ew.rnd$skew, ew.rnd$kurt)
  shm.points = dshimko(r, H, S, K, b,shm.rnd[[1]]$a0, shm.rnd[[1]]$a1, shm.rnd[[1]]$a2)
  bsm.points = dlnorm(K,bsm.rnd$mu, bsm.rnd$zeta)
  
  results.list = list(gb.points, mln.points, ew.points, shm.points, bsm.points, K, 
                      call.strike, iVolInterpol, delta.simple)
  
  return(results.list)
}

