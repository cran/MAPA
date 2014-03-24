# Based on paper:
# Kourentzes N., Petropoulos F. and Trapero J.R. (2014) 
# Improving forecasting by estimating time series structural 
# components across multiple frequencies. 
# International Journal of Forecasting, 30 (2), 291-302
#
# Fotios Petropoulos & Nikolaos Kourentzes (2014)

library(forecast, quietly=TRUE)   # ETS method
library(miscTools, quietly=TRUE)  # Needed for medians by columns
library(parallel, quietly=TRUE)   # Needed for parallel

#-------------------------------------------------
mapa <- function(insample, ppy=NULL, fh=ppy, ifh=1, minimumAL=1, maximumAL=ppy, 
	comb="mean", paral=0, display=0, outplot=1, hybrid=TRUE)
{
# Wrapper to estimate and produce MAPA in- and out-of-sample forecasts
# Uses mapaest and mapafor
# 
# mapa(insample, ppy, fh=ppy, ifh=0, minimumAL=1, maximumAL=ppy, comb="mean", 
#	output="forecast", paral=0, display=0, outplot=1)
#  
# Inputs:
#   insample    = In sample observations of a time series (vector)
#                 If insample == "paper" then it prints paper reference
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If insample is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   fh          = Forecast horizon. Default = ppy
#   ifh         = In-sample forecast horizon. Default = 1
#   minimumAL   = Lower aggregation level to use. Default = 1
#   maximumAL   = Highest aggregation level to use. Default = ppy, maximumAL>1
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1.
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#
# Output:
#   out$infor   = In-sample forecasts
#   out$outfor  = Out-of-sample forecasts
#   out$MSE     = In-sample MSE error
#   out$MAE     = In-sample MAE error
  
  # Paper info
  if (!is.numeric(insample)){
    writeLines("Paper reference: ")
    writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
    writeLines("Improving forecasting by estimating time series structural components")
    writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
                     " 30 (2), 291-302.",sep=""))
    return(invisible())
  }
  
  # Get ppy and fh
  if (is.null(ppy)){
    if (class(insample)=="ts"){
      ppy <- frequency(insample)
      if (is.null(fh)){fh <- ppy}
    } else {
      stop(paste("Input ppy is not given and insample input not ts class.",
                 "Please provide the periods in a season of the time series",
                 "at the sampled frequency."))
    }
  }
  
  # Estimate MAPA
  mapafit <- mapaest(insample, ppy, minimumAL, maximumAL, paral, display)
  # Produce in- and out-of-sample forecasts
  out <- mapafor(insample, mapafit, fh, ifh, comb, outplot, hybrid)
      
}

#-------------------------------------------------
mapasimple <- function(insample, ppy=NULL, fh=ppy, minimumAL=1, maximumAL=ppy, comb="mean", 
	output="forecast", paral=0, display=0, outplot=1, hybrid=TRUE) 
{
# MAPA estimation and forecast
# 
# mapasimple(insample, ppy, fh=ppy, minimumAL=1, maximumAL=ppy, comb="mean", paral=0, 
#	display=0, outplot=1) 
#  
# Inputs:
#   insample    = In sample observations of a time series (vector)
#                 If insample == "paper" then it prints paper reference
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If insample is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   fh          = Forecast horizon. Default = ppy
#   minimumAL   = Lower aggregation level to use. Default = 1, maximumAL>1
#   maximumAL   = Highest aggregation level to use. Default = ppy
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   output      = Type of output. One of "forecast" or "all". Default is "forecast"
#                 If output="all", both forecasts and components estimates per aggregation
#                 level are provided.
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = time series and forecast only;
#                 2 = time series, forecasts and components. Default is 1. 
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE  
#
# Output:
#   forecasts   = Vector with forecasts
#   components  = array with MAPA components

  # Paper info
  if (!is.numeric(insample)){
    writeLines("Paper reference: ")
    writeLines("Kourentzes N., Petropoulos F. and Trapero J.R. (2014)")
    writeLines("Improving forecasting by estimating time series structural components")
    writeLines(paste("across multiple frequencies. International Journal of Forecasting,", 
		                 " 30 (2), 291-302.",sep=""))
    return(invisible())
  }  
  
  # Get ppy and fh
  if (is.null(ppy)){
    if (class(insample)=="ts"){
      ppy <- frequency(insample)
      if (is.null(fh)){fh <- ppy}
    } else {
    stop(paste("Input ppy is not given and insample input not ts class.",
               "Please provide the periods in a season of the time series",
               "at the sampled frequency."))
    }
  }
  
  # Make sure that maximumAL > 1
  if (maximumAL == 1){
    maximumAL <- maximumAL + 1
  }
  
  # Setup parallel processing if required
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", crs, 'cores'))
  }
  
  observations <- length(insample) # number of observations for the in-sample data
  FCs <- array(0, c(maximumAL, 4, fh)) # the forecasts and the forecasts of the components 
									  # will be saved here
   
  # Aggregation and estimation 
  if (paral != 0){  # Parallel run
    FCs_par <- clusterApplyLB(cl, 1:(maximumAL-minimumAL+1), mapasimple.loop, 
      insample=insample, maximumAL=maximumAL, observations=observations, ppy=ppy,
      display=display, fh=fh)  
  } else {          # Serial run
    FCs_par <- vector("list", (maximumAL-minimumAL+1))
    ALvec <- minimumAL:maximumAL
    for (i in 1:(maximumAL-minimumAL+1)){
      FCs_par[[i]] <- mapasimple.loop(ALvec[i], insample, maximumAL, observations, 
        ppy, display, fh)
    }
  }
  
  if (paral == 2){
    # Stop parallel processing
    stopCluster(cl)
  }
  
  # Reshape parallel output
  FCs_par <- do.call(rbind, FCs_par)
  
  FCs <- array(0, c(maximumAL, 4, fh),dimnames=list(paste("AL",minimumAL:maximumAL,sep=""),
    c("ETS","Level","Trend","Season"),paste("t+",1:fh,sep=""))) # the forecasts and the forecasts 
															                              # of the components will be saved here
  for (f in 1:fh){  
    FCs[, , f] <- t(array(FCs_par[,f],c(4,maximumAL)))
  }
  
  # MAPA combination
  combres <- mapacomb(minimumAL,maximumAL,ppy,FCs,comb,observations)
  forecasts <- combres[[1]]
  perm_levels <- combres[[2]]
  perm_seas <- combres[[3]]
  
  # Calculate hybrid model
  if (hybrid==TRUE){
    forecasts <- (FCs[1,1,] + forecasts)/2
  }
  
  # Plot output
  mapaplot(outplot,FCs,maximumAL,perm_levels,perm_seas,observations,insample,forecasts,fh,comb)
  
  # Construct output
  if (output=="forecast"){
    forecasts
  } else {
    list(forecast=forecasts,components=FCs) 
  }
  
}

#-------------------------------------------------
mapafor <- function(insample, mapafit, fh=-1, ifh=1, comb="mean", outplot=1, hybrid=TRUE) {
# MAPA in- and out-of-sample forecast
# 
# mapafor(insample, mapafit, fh=ppy, ifh=0, comb="mean", outplot=1) 
#  
# Inputs:
#   insample    = In sample observations of a time series (vector)
#   mapafit     = Fitted MAPA model (from mapaest)
#   fh          = Forecast horizon. Default = ppy
#   ifh         = In-sample forecast horizon. Default = 1
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1. 
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#
# Output:
#   out$infor   = In-sample forecasts
#   out$outfor  = Out-of-sample forecasts
#   out$MSE     = In-sample MSE error
#   out$MAE     = In-sample MAE error
  
  observations <- length(insample) # number of observations for the in-sample data
  
  ppy <- as.numeric(mapafit[1,12])
  if (fh == -1){
    fh <- ppy
  }
  
  # In-sample MAPA
  if (ifh>0){
    infor <- array(NA,c(ifh,observations),dimnames=list(paste("t+",1:ifh,sep="")))
    for (i in (ppy):(observations-1)){
      inobs <- as.matrix(insample[1:i])
      infor[, i+1] <- mapacalc(inobs, mapafit, fh=ifh, comb, output="forecast", outplot=0, hybrid) 
      # Crop out-of-sample predictions
      if ((i+ifh)>observations){
        k <- (i+ifh) - observations
        infor[(ifh-k+1):ifh, i+1] <- rep(NA,k)
      }
    }    
  } else {
    infor <- NULL
  }
  
  # Out-of-sample MAPA
  if (fh>0){
    outfor <- mapacalc(insample, mapafit, fh, comb, output="forecast", outplot=0, hybrid) 
  } else {
    outfor <- NULL
  }
  
  # Produce plot
  if (outplot==1){
    # Find min max
    if (is.null(outfor)){
      ymax <- max(insample)
      ymin <- min(insample)
      ymax <- ymax + 0.1*(ymax-ymin)
      ymin <- ymin - 0.1*(ymax-ymin)      
    } else {
      ymax <- max(c(max(outfor),max(insample)))
      ymin <- min(c(min(outfor),min(insample)))
      ymax <- ymax + 0.1*(ymax-ymin)
      ymin <- ymin - 0.1*(ymax-ymin)
    }
    plot(1:observations,insample,type="l",col="blue", xlab="", ylab="", main="Forecast", 
		xlim <- c(1, observations+fh), ylim=c(ymin,ymax))
    # In-sample
    if (ifh>0){
      if (ifh==1){
        lines(infor[1,],col="red")
      } else {
        # clrs = rainbow(observations-ppy)
        for (i in (ppy):(observations-1)){
          lines((i):(i+ifh-1),infor[,i],col="red")
        }
      }
    }
    # Out-of-sample
    if (ifh == 0){
      lines(observations:(fh+observations),c(insample[observations],outfor),col="red")
    } else if (ifh == 1){
      lines(observations:(fh+observations),c(infor[1,observations],outfor),col="red")
    } else {
      lines((observations+1):(fh+observations),outfor,col="red")
    }
  }
  
  # Calculate in-sample errors
  if (ifh == 1) {
    resid <- insample - t(infor)
    MSE <- mean(resid^2, na.rm=TRUE)
    MAE <- mean(abs(resid), na.rm=TRUE)
  } else if (ifh > 1) {
    MSE <- array(0,c(ifh,1),dimnames=list(paste("t+",1:ifh,sep=""),"MSE"))
    MAE <- array(0,c(ifh,1),dimnames=list(paste("t+",1:ifh,sep=""),"MAE"))
    for (h in 1:ifh) {
      resid <- insample[h:observations] - infor[h, 1:(observations-h+1)]
      MSE[h] <- mean(resid^2, na.rm=TRUE)
      MAE[h] <- mean(abs(resid), na.rm=TRUE)
    }
  } else {
    MSE <- NULL
    MAE <- NULL
  }
  
  # Construct output
  output <- list(infor=infor,outfor=outfor,MSE=MSE,MAE=MAE)
  
}

#-------------------------------------------------
mapaest <- function(insample, ppy=NULL, minimumAL=1, maximumAL=ppy, paral=0, display=0, outplot=1) {
# Estimate MAPA for a time series  
# 
# mapaest(insample, ppy, minimumAL=1, maximumAL=ppy, paral=0, display=0) 
#  
# Inputs:
#   insample    = In sample observations of a time series (vector)
#   ppy         = Periods in a season of the time series at the sampled frequency.
#                 If insample is a ts object then this is taken from its frequency,
#                 unless overriden. 
#   minimumAL   = Lower aggregation level to use. Default = 1, maximumAL>1
#   maximumAL   = Highest aggregation level to use. Default = ppy
#   paral       = Use parallel processing. 0 = no; 1 = yes (requires initialised cluster); 
#                 2 = yes and initialise cluster. Default is 0.
#   display     = Display calculation progress in console. 0 = no; 1 = yes. Default is 0.
#   outplot     = Provide output plot. 0 = no; 1 = yes. Default is 1.   
#
# Output:
#   mapafit     = Estimated MAPA model structure

  # Get ppy and maximumAL
  if (is.null(ppy)){
    if (class(insample)=="ts"){
      ppy <- frequency(insample)
      if (is.null(maximumAL)){maximumAL <- ppy}
    } else {
      stop(paste("Input ppy is not given and insample input not ts class.",
                 "Please provide the periods in a season of the time series",
                 "at the sampled frequency."))
    }
  }  
  
  # Make sure that maximumAL > 1
  if (maximumAL == 1){
    maximumAL = maximumAL + 1
  }
  
  # Setup parallel processing if required
  if (paral == 2){
    crs <- detectCores()
    cl <- makeCluster(getOption("cl.cores", crs))
    writeLines(paste("Running with", crs, 'cores'))
  }
  
  observations <- length(insample) # number of observations for the in-sample data
  
  # Aggregation and estimation 
  if (paral != 0){  # Parallel run
    mapafit <- clusterApplyLB(cl, 1:(maximumAL-minimumAL+1), mapaest.loop, 
      insample=insample, maximumAL=maximumAL, observations=observations, ppy=ppy,
      display=display)  
  } else {          # Serial run
    mapafit <- vector("list", (maximumAL-minimumAL+1))
    ALvec <- minimumAL:maximumAL
    for (i in 1:(maximumAL-minimumAL+1)){
      mapafit[[i]] <- mapaest.loop(ALvec[i], insample, maximumAL, observations, 
        ppy, display)
    }
  }
    
  if (paral == 2){
    # Stop parallel processing
    stopCluster(cl)
  }

  # Process output
  mapafit <- do.call(rbind, mapafit) # Necessary for clusterApplyLB function
  rownames(mapafit) <- paste("AL",minimumAL:maximumAL,sep="")

  # mapafit <- mapafit[, c(20,19,11,12,13,14)]
  
  # Plot model selection summary
  if (outplot == 1){
    comps <- array(0,c(maximumAL,5))
    for (AL in minimumAL:maximumAL){
      components <- mapafit[[AL, 14]]
      # Error term
      if (components[1]=="A"){
        comps[AL,1] <- 1
      } else {
        comps[AL,1] <- 2
      }
      # Trend term
      if (components[2]=="A"){
        comps[AL,2] <- 1
      } else {if (components[2]=="M"){
        comps[AL,2] <- 2
      } else
        comps[AL,2] <- 0
      }
      # Season term
      if (components[3]=="A"){
        comps[AL,3] <- 1
      } else {if (components[3]=="M"){
        comps[AL,3] <- 2
      } else
        comps[AL,3] <- 0
      }
      # Damped tem
      if (components[4]==TRUE){
        comps[AL,4] <- 1
      }
      comps[AL,5] <- mapafit[[AL,20]]
    }
    # Use only appropriate ones
    comps <- comps[unlist(mapafit[, 19])==TRUE, ]
    comps[, 2] <- comps[, 2] + 0.5*comps[, 4]
    image(min(comps[,5]):max(comps[,5]), 1:3, comps[,1:3], axes=FALSE, col=rev(heat.colors(5)), 
		  ylab="Components", xlab="Aggregation Level", main="ETS components")
    axis(2, at=1:3, labels=list("Error","Trend","Season"))
    axis(1, at=min(comps[,5]):max(comps[,5]))
    
    for (i in 1:4){
      for (AL in min(comps[,5]):(max(comps[,5])+1)){
        if (i==1){
          lines(c(AL-0.5,AL-0.5),c(0,4),col="black")
        }
        if (i<4 & AL<=max(comps[,5])){
          if (i==2 & comps[AL,4]==TRUE){
            damp <- "d"
          } else {
            damp <- NULL
          }
          text(AL,i,paste(mapafit[[AL,14]][i],damp,sep=""))
        }
      }
      lines(c(min(comps[,5])-0.5,max(comps[,5])+0.5),c(i-0.5,i-0.5),col="black")
    }
  }
  
  # Return output
  mapafit
  
}

#-------------------------------------------------
mapacalc <- function(insample, mapafit, fh=0, comb="mean", output="forecast", 
	outplot=0, hybrid=TRUE) 
{
# Calculation of MAPA forecasts
# 
# mapacalc(insample, mapafit, fh=0, comb="mean", output="forecast", outplot=0) 
#  
# Inputs:
#   insample    = In sample observations of a time series (vector)
#   mapafit     = Fitted MAPA model (from mapaest)
#   fh          = Forecast horizon. Default = ppy
#   comb        = Combination operator. One of "mean" or "median". Default is "mean"
#   output      = Type of output. One of "forecast" or "all". Default is "forecast"
#                 If output="all", both forecasts and components estimates per aggregation
#                 level are provided.
#   outplot     = Provide output plot. 0 = no; 1 = time series and forecast only;
#                 2 = time series, forecasts and components. Default is 1. 
#   hybrid      = Provide hybrid forecasts, as in Kourentzes et al. paper. Default is TRUE
#
# Output:
#   forecasts   = Vector with forecasts
#   components  = array with MAPA components
  
  # Get settings from mapafit
  ALs <- as.numeric(mapafit[mapafit[,19]==TRUE, 20])
  minimumAL <- min(ALs)
  maximumAL <- max(ALs)
  ppy <- as.numeric(mapafit[1,12])
  if (fh==0){
    fh <- ppy
  } 
  
  observations <- length(insample) # number of observations for the in-sample data
  
  FCs <- array(0, c(maximumAL, 4, fh),dimnames=list(paste("AL",minimumAL:maximumAL,sep=""),
	c("ETS","Level","Trend","Season"),paste("t+",1:fh,sep=""))) # the forecasted components 
																# are saved here
  
  # MAPA forecast
  for (AL in minimumAL:maximumAL){
    
    q <- observations %/% AL # observation in the aggregated level
    r <- observations %% AL  # observation to discard from the beginning of the series
    ppyA <- ppy %/% AL       # periods per year for the aggregated level
    if (ppy %% AL != 0){
      ppyA <- 1
    }
    
    # Aggregation
    insampleA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    for (j in 1:q){                 # calculate the aggregate values
      insampleA[j] <- mean(insample[(r+1+(j-1)*AL):(r+j*AL)])
    }
    
    ats <- ts(insampleA, frequency = ppyA) # create the time series

# This part of code used to call forecast:::pegelsresid.C
# This is now obsolete with the new forecast package that has use.initial.values=TRUE
#
#     # Extarct ets fit from mapafit
#     components <- mapafit[[AL,14]]
#     errortype <- components[1]
#     trendtype <- components[2]
#     seasontype <- components[3]
#     damped <- as.logical(components[4])
#     param <- mapafit[[AL,11]]
#     pnames <- names(param)
#     alpha <- param[pnames=="alpha"]
#     if (sum(pnames=="beta")>0) {
#       beta <- param[pnames=="beta"]
#     } else {
#       beta <- NULL
#     }
#     if (sum(pnames=="gamma")>0) {
#       gamma <- param[pnames=="gamma"]
#     } else {
#       gamma <- NULL
#     }
#     if (sum(pnames=="phi")>0) {
#       phi <- param[pnames=="phi"]
#     } else {
#       phi <- NULL
#     }
#     
#     # Prepare initial states
#     np <- length(param)
#     nstate <- np - length(c(alpha,beta,gamma,phi))
#     initstate <- param[(np-nstate+1):np]
#     # Add extra state as per ets.R
#     if(seasontype!="N")
#       initstate <- c(initstate, ppyA*(seasontype=="M")-sum(initstate[(2+(trendtype!="N")):nstate]))
#     
#     # Calculate ets states using mapafit results
#     ats.fit <- forecast:::pegelsresid.C(ats, ppyA, initstate, errortype, trendtype,
# 		  seasontype, damped, alpha, beta, gamma, phi)
#     ats.fit$components <- components
    
    # ETS based calculation
    param <- mapafit[[AL,11]]
    pnames <- names(param)
    if (sum(pnames=="phi")>0) {
      phi <- param[pnames=="phi"]
    } else {
      phi <- NULL
    }
    
    AL.fit <- structure(mapafit[AL,1:18],class="ets")
    ats.fit <- ets(ats, AL.fit, use.initial.values=TRUE)
    
    # Transalte ets states for MAPA
    FCs_temp <- statetranslate(ats.fit,AL,fh,q,ppyA,phi,1)
    
    # Return MAPA components
    FCs[AL, , ] <- FCs_temp
    
  }
  
  # MAPA combination
  combres <- mapacomb(minimumAL,maximumAL,ppy,FCs,comb,observations)
  forecasts <- combres$forecasts
  perm_levels <- combres$perm_levels
  perm_seas <- combres$perm_seas

  # Calculate hybrid model
  if (hybrid==TRUE){
    forecasts <- (FCs[1,1,] + forecasts)/2
  }  
  
  # Plot output
  mapaplot(outplot,FCs,maximumAL,perm_levels,perm_seas,observations,insample,forecasts,fh,comb)

  # Construct output
  if (output=="forecast"){
    forecasts
  } else {
    list(forecast=forecasts,components=FCs) 
  }
  
}

#-------------------------------------------------
statetranslate <- function(fit,AL,fh,q,ppyA,phi,fittype){
# This function prepares ets states for MAPA combination
# It extrapolates from last states the forecasts and translates to additive
  
  FCs_temp <- array(0, c(4, fh))
  
  fhA <- (fh %/% AL) + 1   # forecast horizon for the aggregated level
  
  # Estimates for the Level Component
  FCs_temp[2, ] <- as.numeric(rep(rep(fit$states[q+1, 1], fhA), each=AL)[1:fh])
  
  # Estimates for the Trend Component
  if (fit$components[2]=="N"){ # no trend
    FCs_temp[3, ] <- 0
    b = 0 # indicates that there is no trend
  } else if (fit$components[2]=="A"){ # additive trend
    if (fit$components[4]=="FALSE"){ 
      FCs_temp[3, ] <- as.numeric(rep(fit$states[q+1, 2] * (1:fhA), each=AL))[1:fh]
    } else { # additive damped trend
      # We divide with phi because of an internal calculation for 
	  # the damped trend in the ETS package
      FCs_temp[3, ] <- as.numeric(rep(cumsum((fit$states[q+1, 2]/phi)* phi^(1:fhA)), each=AL))[1:fh]
    }
    b <- 1 # indicates that there is trend
  } else {
    if (fit$components[4]=="FALSE"){ # multiplicative trend
      FCs_temp[3, ] <- as.numeric(rep((fit$states[q+1,2]^(1:fhA)-1), each=AL)[1:fh] * FCs_temp[2,])
    } else { # multiplicative damped trend
      # We divide with phi because of an internal calculation for the damped trend in the ETS package
      FCs_temp[3, ] <- as.numeric(rep((((fit$states[q+1,2] ^ (1/phi)) ^ cumsum(phi^(1:fhA)))-1), 
		each=AL)[1:fh] * FCs_temp[2, ])
    }
    b <- 1 # indicates that there is trend
  }
  
  # Estimates for the Seasonal Component 
  if (fit$components[3]=="N"){ # no seasonality
    FCs_temp[4, ] <- 0
  } else if (fit$components[3]=="A"){ # additive seasonality
    FCs_temp[4, ] <- as.numeric(rep(rep(rev(fit$states[q+1,(2+b):(ppyA+1+b)]), fhA), each=AL))[1:fh]
  } else { # multiplicative seasonality
    FCs_temp[4, ] <- as.numeric((rep(rep(rev(fit$states[q+1,(2+b):(ppyA+1+b)]), fhA),
		each=AL)[1:fh] - 1)) * (FCs_temp[2, ] + FCs_temp[3, ])
  }
  
  # fittype identifies if information is comming from ets or mapafit
  if (fittype==1){
    # Recreate ETS forecasts
    if (fh != 1) {
      FCs_temp[1, ] <- colSums(FCs_temp[2:4,])     
    } else {
      FCs_temp[1, ] <- sum(FCs_temp[2:4,])   
    }
  }

  # Return output
  FCs_temp
  
}

#-------------------------------------------------
mapacomb <- function(minimumAL,maximumAL,ppy,FCs,comb,observations){
# This function combines the translated ets states
  
  perm_levels <- array(0, maximumAL) # permitted levels due to ETS implementation (observations>=4)
  perm_seas <- array(0, maximumAL)   # permitted seasonalities
  for (AL in minimumAL:maximumAL){
    if (observations %/% AL >=4){
      perm_levels[AL] <- 1
    }
    if ((ppy %% AL == 0) & (AL<ppy)){
      perm_seas[AL] <- 1
    }
  }
  
  if (dim(FCs)[3] != 1){ # Forecast multiple steps ahead
    if (comb=="mean"){ # alternative averaging operators
      forecasts <- colSums(rbind(colMeans(FCs[perm_levels==1, 2, ]),
		    colMeans(FCs[perm_levels==1, 3, ]),colMeans(FCs[(perm_levels==1 & perm_seas==1), 4, ])),
		    na.rm=TRUE) # MAPA(mean) forecasts
    } else {
      forecasts <- colSums(rbind(colMedians(FCs[perm_levels==1, 2, ]),
		    colMedians(FCs[perm_levels==1, 3,]),colMedians(FCs[(perm_levels==1 & perm_seas==1), 4, ])),
		    na.rm=TRUE) # MAPA(median) forecasts
    }
  } else {
    if (comb=="mean"){ # alternative averaging operators
      forecasts <- colSums(rbind(mean(FCs[perm_levels==1, 2, ]),mean(FCs[perm_levels==1, 3, ]),
		    mean(FCs[(perm_levels==1 & perm_seas==1), 4, ])), na.rm=TRUE) # MAPA(mean) forecasts
    } else {
      forecasts <- colSums(rbind(median(FCs[perm_levels==1, 2, ]),median(FCs[perm_levels==1, 3, ]),
		    median(FCs[(perm_levels==1 & perm_seas==1), 4, ])), na.rm=TRUE) # MAPA(mean) forecasts
    }
  }
  
  # Return output
  list(forecasts=forecasts,perm_levels=perm_levels,perm_seas=perm_seas)
}

#-------------------------------------------------
mapaplot <- function(outplot,FCs,maximumAL,perm_levels,perm_seas,observations,
	insample,forecasts,fh,comb)
{
# Produce MAPA forecast & components plot 
# outplot == 0, no plot; == 1 series plot; == 2 component plot
  
  # Make sure that maximumAL > 1
  if (maximumAL == 1){
    maximumAL <- maximumAL + 1
  }
  
  if (outplot > 0){
    FClevel <- FCs[perm_levels==1, 2, ]
    FCtrend <- FCs[perm_levels==1, 3, ]
    if (sum(perm_levels==1 & perm_seas==1)!=0){
      FCseason <- FCs[(perm_levels==1 & perm_seas==1), 4, ]
    } else {
      FCseason <- NULL
    }
    clrs <- rainbow(maximumAL)
    if (outplot == 2){
      if (!is.null(FCseason)){
        layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
      } else {
        layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
      }
    } else {
      layout(matrix(1, 1, 1, byrow = TRUE))
    }
    # Find min max
    ymax <- max(c(max(forecasts),max(insample)))
    ymin <- min(c(min(forecasts),min(insample)))
    ymax <- ymax + 0.1*(ymax-ymin)
    ymin <- ymin - 0.1*(ymax-ymin)
    # Plot prediction
    plot(1:observations,insample, type="l", col="blue", xlab="", ylab="", 
		  main="Forecast", xlim=c(1, observations+fh), ylim=c(ymin,ymax))
    lines((observations):(observations+fh),c(insample[observations],forecasts), col="red")
    if (outplot == 2){
      # Plot Level
      plot(FClevel[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Level")
      for (i in 2:sum(perm_levels)){
        lines(FClevel[i, ], type="l", col=clrs[i])
      }
      if (comb=="mean"){
        lines(colMeans(FClevel), type="l", col="black", lwd=2)
      } else {
        lines(colMedians(FClevel), type="l", col="black", lwd=2)
      }
      # Plot trend
      plot(FCtrend[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Trend", 
		    ylim=c(min(FCtrend),max(FCtrend)))
      for (i in 2:sum(perm_levels)){
        lines(FCtrend[i, ], type="l", col=clrs[i])
      }
      if (comb=="mean"){
        lines(colMeans(FCtrend), type="l", col="black", lwd=2)
      } else {
        lines(colMedians(FCtrend), type="l", col="black", lwd=2)
      }
      # Plot season
      if (!is.null(FCseason)){
        plot(FCseason[1, ], type="l", col=clrs[1], xlab="", ylab="", main="Season")
        for (i in 2:sum(perm_levels==1 & perm_seas==1)){
          lines(FCseason[i, ], type="l", col=clrs[i])
        }
        if (comb=="mean"){
          lines(colMeans(FCseason), type="l", col="black", lwd=2)
        } else {
          lines(colMedians(FCseason), type="l", col="black", lwd=2)
        }
      }
    }
  }
}

#-------------------------------------------------
mapasimple.loop <- function(AL, insample, maximumAL, observations, ppy, display,fh){
  # Internal function for mapasimple estimation and forecast iterations
  
  # Console output
  if (display==1){
    txtc <- paste("Aggregation level: ",AL,"/",maximumAL,
                  " (",round(100*AL/maximumAL,2),"%)",sep="")
    cat(txtc)
  }
  
  FCs_temp <- array(0, c(4, fh))
  
  q <- observations %/% AL # observation in the aggregated level
  r <- observations %% AL  # observation to discard from the beginning of the series
  fhA <- (fh %/% AL) + 1   # forecast horizon for the aggregated level
  ppyA <- ppy %/% AL       # periods per year for the aggregated level
  if (ppy %% AL != 0){
    ppyA <- 1
  }
  
  if (q >= 4){
    # Aggregation
    insampleA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    for (j in 1:q){                # calculate the aggregate values
      insampleA[j] <- mean(insample[(r+1+(j-1)*AL):(r+j*AL)])
    }
    
    ats <- ts(insampleA, frequency = ppyA) # create the time series
    
    # Fit ETS
    ats.fit <- ets(ats)
    ats.fcast <- forecast.ets(ats.fit,h=fhA)
    
    # Translate ets states for MAPA
    phi <- ats.fit$par[names(ats.fit$par)=="phi"]
    FCs_temp <- statetranslate(ats.fit,AL,fh,q,ppyA,phi,0)
    
    # ets forecast on original frequency
    FCs_temp[1, ] <- rep(ats.fcast$mean[1:fhA], each=AL)[1:fh]
  }
  
  # Update console display
  if (display==1){
    nc <- nchar(txtc)
    cat(rep("\r",nc))
    cat(rep(" ",nc))
    cat(rep("\r",nc))
  }
  
  # Return FCs_temp for foreach
  FCs_temp
}

#-------------------------------------------------
mapaest.loop <- function(AL, insample, maximumAL, observations, ppy, display){ 
  # Internal function for running a single loop in mapaest
  
  # Console output
  if (display==1){
    txtc <- paste("Aggregation level: ",AL,"/",maximumAL,
                  " (",round(100*AL/maximumAL,2),"%)",sep="")
    cat(txtc)
  }
  
  q <- observations %/% AL # observation in the aggregated level
  r <- observations %% AL  # observation to discard from the beginning of the series
  ppyA <- ppy %/% AL       # periods per year for the aggregated level
  if (ppy %% AL != 0){
    ppyA <- 1
  }
  
  if (q >= 4){
    # Aggregation
    insampleA <- array(0, dim=c(q)) # in-sample aggregated values will be saved here
    for (j in 1:q){                 # calculate the aggregate values
      insampleA[j] <- mean(insample[(r+1+(j-1)*AL):(r+j*AL)])
    }
    
    ats <- ts(insampleA, frequency = ppyA) # create the time series
    
    # Fit ETS
    fit <- ets(ats)
    fit$use <- TRUE
  } else {
    fit <- NULL
    fit$use <- FALSE
  }
  
  fit$AL <- AL
  
  # Update console display
  if (display==1){
    nc = nchar(txtc)
    cat(rep("\r",nc))
    cat(rep(" ",nc))
    cat(rep("\r",nc))
  }
  
  # Return loop result
  rbind(fit)
}