#' Track filtering and simulation
#'
#' Provides facilities for filtering and simulating tracks from
#' telemetry to support the Marseille RAATD workshop.
#'
#' @name Marseille-package
#' @docType package
#' @author Ben Raymond, Ian Jonsen, Ryan Reisinger and Simon Wotherspoon
NULL


##' Correlated Random Walk Filter
##'
##' Fit a correlated random walk to filter a track and predict
##' locations on a regular time step.
##'
##' The input track is given as a dataframe where each row is an
##' observed location, with columns
##' \describe{
##' \item{'date'}{observation time (as GMT POSIXct),}
##' \item{'lon'}{observed longitude,}
##' \item{'lat'}{observed latitude,}
##' \item{'lc'}{location class.}
##' }
##'
##' The TMB parameter list can be specified directly with the
##' \code{parameters} argument. Otherwise suitable parameter values
##' are estimated from predicted locations estimated by fitting loess
##' smooths to the raw observations.
##'
##' If \code{extrap} is \code{TRUE}, the last observations occur after
##' the last predicted location and the last fitted locations are
##' extrapolations, otherwise the final observations occur before the
##' final predicted locations and all fitted locations are
##' interpolated.
##'
##' The filtering model assumes the errors in longitude and latitude
##' are proportional to scale factors determined by the location
##' class. The scale factors are speficied through the \code{aes}
##' argument. By default the function uses the same scaling factors
##' for location accuracy as used in \pkg{crawl} for ARGOS data.
##'
##' @title Correlated Random Walk Filter
##' @param data A dataframe representing the track (see details)
##' @param subset Logical vector indicating which rows of the data
##'   frame be kept in the filtering
##' @param tstep The time step to predict to (in days)
##' @param nu The degrees of freedom parameter for the location errors
##' @param gamma The autocorrelation parameter used to estimate
##'   initial parameters.
##' @param span The span parameter for the loess fits used to estimate
##'   initial locations.
##' @param verbose Enable tracing information
##' @param optim The function used to minimize the negative log
##'   likelihood.
##' @param extrap If TRUE, the final predicted state occurs
##'   immediately before the last observation, otherwise the final
##'   predicted state occurs immediately after the last observation.
##' @param parameters The TMB parameter list.
##' @param esf The error scale factors for the location classes.
##' @param common.tau Should a common tau parameter be fitted for
##'   longitude and latitude.
##' @return a list with components
##'   \item{\code{predicted}}{a dataframe of predicted locations}
##'   \item{\code{fitted}}{a dataframe of fitted locations}
##'   \item{\code{par}}{model parameter summary}
##'   \item{\code{data}}{the input dataframe}
##'   \item{\code{subset}}{the input subset vector}
##'   \item{\code{tstep}}{the prediction time step}
##'   \item{\code{common.tau}}{has a common tau been fitted for lon and lat}
##'   \item{\code{opt}}{the object returned by the optimizer}
##'   \item{\code{tmb}}{the \pkg{TMB} object}
##' @useDynLib Marseille
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats sd cov loess loess.control predict nlminb optim
##' @export
dcrw <- function(data,subset=rep(TRUE,nrow(data)),
                 tstep=6/24, nu=10, gamma=0.5, span=0.1, verbose=FALSE,
                 optim=c("nlminb","optim"), extrap=FALSE, parameters=NULL,
                 esf=aesfCRAWL(),common.tau=FALSE) {

  optim <- match.arg(optim)
  d <- data[subset,]

  ## TMB - create data list
  d$date <- as.POSIXct(d$date, "%Y-%m-%d %H:%M:%S", tz="GMT")
  d$lc <- factor(d$lc,levels=levels(esf$lc),ordered=TRUE)
  d$lon <- unwrapLon(d$lon)
  ## Merge error multiplication factors
  d <- merge(d, esf, by="lc", all.x=TRUE)
  d <- d[order(d$date),]

  ## Interpolation indices and weights
  dt <- tstep*86400
  tms <- (as.numeric(d$date)-as.numeric(d$date[1]))/dt
  index <- floor(tms)
  if(extrap) index <- pmax(index,max(index)-1)
  weights <- 1-(tms-index)

  tmbdata <- list(y=cbind(d$lon,d$lat),idx=index,w=weights,
                  K=cbind(d$lonESF,d$latESF),nu=nu)

  ## Predict track from loess smooths
  fit.lon <- loess(lon ~ as.numeric(date), data=d, span=span,
                     na.action="na.exclude",control = loess.control(surface = "direct"))
  fit.lat <- loess(lat ~ as.numeric(date), data=d, span=span,
                     na.action="na.exclude",control = loess.control(surface = "direct"))
  ## Predict track, increments and stochastic innovations
  ts <- seq(d$date[1],by=dt,length.out=max(index)+2)
  xs <- cbind(predict(fit.lon,newdata=data.frame(date=as.numeric(ts))),
              predict(fit.lat,newdata=data.frame(date=as.numeric(ts))))


  ## TMB - create parameter list
  if(is.null(parameters)) {
    ## Estimate increments and stochastic innovations
    ds <- xs[-1,] - xs[-nrow(xs),]
    es <- ds[-1,] - gamma*ds[-nrow(ds),]

    ## Estimate components of variance
    V <- cov(es)
    sigma <- sqrt(diag(V))
    rho <- V[1,2]/prod(sqrt(diag(V)))
    if(common.tau)
      tau <- sd(c(fit.lon$residuals/d$lonESF,fit.lat$residuals/d$latESF))
    else
      tau <- c(sd(fit.lon$residuals/d$lonESF),sd(fit.lat$residuals/d$latESF))

    parameters <- list(theta=0,
                       l_gamma=log(gamma/(1-gamma)),
                       l_sigma=log(pmax(1.0E-8,sigma)),
                       l_rho=log((1+rho)/(1-rho)),
                       l_tau=log(pmax(1.0E-8,tau)),
                       x=xs)
  }

  ## TMB - create objective function
  obj <- MakeADFun(tmbdata,parameters,random="x",DLL="Marseille",silent=!verbose)
  obj$env$inner.control$trace <- verbose
  obj$env$tracemgc <- verbose

  ## Minimize objective function
  opt <- suppressWarnings(
    switch(optim,
           nlminb=nlminb(obj$par, obj$fn, obj$gr),
           optim=do.call("optim", obj)))

  ## Parameters, states and the fitted values
  rep <- sdreport(obj)
  fxd <- summary(rep,"report")
  rdm <- matrix(summary(rep,"random"),length(ts),4,
                dimnames=list(NULL,c("lon","lat","lon.se","lat.se")))
  ftd <- weights*rdm[index+1,1:2]+(1-weights)*rdm[index+2,1:2]

  list(predicted=cbind(date=ts,as.data.frame(rdm)),
       fitted=cbind(date=d$date,as.data.frame(ftd)),
       par=fxd,data=data,subset=subset,tstep=tstep,
       common.tau=common.tau,opt=opt,tmb=obj)
}



##' Wrappng Locations Around the Dateline
##'
##' These functions wrap and unwrap a sequence of longitudes around
##' the dateline.
##'
##' The \code{wrapLon} function wraps the longitudes back into the
##' interval [lmin,lmin+360).  The \code{unwrapLon} function unwraps a
##' sequence of longitudes so the the initial point lies in
##' [lmin,lmin+360), but the subsequent longitudes in the sequence may
##' wander outside that range.
##'
##' @title Dateline adjustment
##' @param lon a vector of longitudes
##' @param lmin western boundary for wrapped longitudes
##' @return a vector of longitudes
##' @export
wrapLon <- function(lon,lmin=-180)
  (lon-lmin)%%360+lmin

##' @rdname wrapLon
##' @export
unwrapLon <- function(lon,lmin=-180)
  cumsum(c(wrapLon(lon[1],lmin),wrapLon(diff(lon))))



##' ARGOS Error Scale Factors for Location Classes
##'
##' These are the multiplication factors used to scale location
##' accuracy for each location class that are used in \pkg{crawl}.
##' @title ARGOS Error Scale Factors
##' @return A dataframe with columns
##' \item{\code{lc}}{ARGOS location class}
##' \item{\code{lonESF}}{error scaling factor for longitude}
##' \item{\code{latESF}}{error scaling factor for latitude}
##' @export
aesfCRAWL <- function() {
  data.frame(lc=factor(c("3", "2", "1", "0", "A", "B"),
                       levels=c("3", "2", "1", "0", "A", "B"),ordered=TRUE),
             lonESF=c(1,1.54,3.72,23.90,13.51,44.22),
             latESF=c(1,1.29,2.55,103.70,14.99,32.53))
}


##' Generate Transition and Covariance Matrices for a DCRW model
##'
##' Generates the transition and covariance matrices corresponding to
##' the dcrw state space model. Currently no polar adjustment is made.
##'
##' @title DCRW Movement Model
##' @param fit a fitted \code{dcrw} object
##' @return A list containing
##' \item{\code{A}}{the transition matrix}
##' \item{\code{Q}}{the covariance matrix}
##' @export
dcrwModel <- function(fit) {

  ## Transition matrix A = [I, I; 0 gT]
  theta <- fit$par[1,1]
  gamma <- fit$par[2,1]
  A <- matrix(0,4,4)
  A[1:2,1:2] <- diag(1,2,2)
  A[1:2,3:4] <- diag(1,2,2)
  A[3:4,3:4] <- matrix(gamma*c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)

  ## Innovation variance Q = [0, 0; 0, Sigma]
  sigma <- fit$par[3:4,1]
  rho <- fit$par[5,1]
  Q <- matrix(0,4,4)
  Q[3:4,3:4] <- diag(sigma,2,2)%*%matrix(c(1,rho,rho,1),2,2)%*%diag(sigma,2,2)

  list(A=A,Q=Q)
}

##' Generate new tracks from a fitted DCRW model
##'
##' Given a fitted dcrw model this function generates a new track of
##' the same length that coincides with the fitted track at the start
##' point and optionally other specified points along the estimated
##' track.
##'
##' Locations from the fitted track can be marked as fixed with the
##' \code{fixed} argument. In the current implementation the first
##' location must always be fixed.  The \code{fixed.err} parameter
##' specifies the covariance of the error in the fixed points. If this
##' parameter is \code{NULL}, the covariance defaults to the
##' innovation covariance \code{Q} returned by \code{dcrwModel}.
##'
##' Additional constraints can be placed on the path by rejection
##' sampling through the function \code{point.check}.  This function
##' must accept a time, longitude and latitude and return a boolean
##' indicating whether the point is acceptable.  For example, the
##' track can be constrained to the ocean by supplying a
##' \code{point.check} function that compares the state to a land mask
##' and returns \code{FALSE} for locations on land.
##'
##' Tracks are simulated in the plane.  There is is no polar
##' correction and without a suitable \code{point.check} function,
##' there is nothing prevent the track extending outside the [-90,90]
##' latitudinal limits.
##'
##' @title Regressive bridge sampler
##' @param fit an object returned by \code{dcrw}
##' @param fixed a logical vector indicating which locations in the
##'   template path are to be held fixed.
##' @param fixed.err The covariance matrix for fixed points.
##' @param point.check function that accepts a time, longitude and
##'   latitude and returns boolean indicating whether the state is
##'   acceptable.
##' @return a dataframe representing the simulated track
##'   \item{\code{date}}{prediction times as POSIXct.}
##'   \item{\code{lon}}{predicted longitude}
##'   \item{\code{lat}}{predicted latitude}
##' @importFrom stats rnorm
##' @export
dcrwSimulate <- function(fit, fixed=rep(c(TRUE,FALSE,TRUE),c(1,nrow(fit$predicted)-2,1)),
                         fixed.err=NULL, point.check=function(tm,lon,lat) TRUE) {

  ## Extract model matrices
  model <- dcrwModel(fit)
  A <- model$A
  Q <- model$Q+diag(1.0E-8,4,4)
  if(is.null(fixed.err)) fixed.err <- Q

  ## First state must be fixed
  fixed[1] <- TRUE

  ## Times and matrix of states
  ts <- fit$predicted$date
  xs <- cbind(unname(as.matrix(fit$predicted[,c("lon","lat")])),0,0)
  n <- nrow(xs)

  ## Prior mean and variance for each state
  ms <- xs
  Vs <- array(0,c(nrow(Q),ncol(Q),n))

  ## Forward pass - generate priors from movement model
  for(k in 1:n)
    if(fixed[k]) {
      ms[k,3:4] <- 0
      Vs[,,k] <- fixed.err
    } else {
      ms[k,] <- A%*%ms[k-1,]
      Vs[,,k] <- A%*%Vs[,,k-1]%*%t(A)+Q
    }

  ## Reverse pass - recursively sample with a Kalman/Regression step
  ## starting from x[k0,]
  sample <- function(k0) {
    x <- ms[k0,] + drop(rnorm(4)%*%chol(Vs[,,k0]))
    xs[k0,] <<- x
    for(k in (k0-1):1) {
      ## Kalman gain
      K <- Vs[,,k]%*%t(A)%*%solve(A%*%Vs[,,k]%*%t(A)+Q)
      ## Mean, variance update
      mu <- ms[k,] + drop(K%*%(x-A%*%ms[k,]))
      V <- Vs[,,k] - K%*%A%*%Vs[,,k]
      ##W <- (diag(1,4,4)-K%*%A)
      ##V <- tcrossprod(W%*%Vs[,,k],W)+tcrossprod(K%*%Q,K)
      R <- chol(V)
      ## point.check/rejection loop
      for(r in 1:100) {
        x <- mu + drop(rnorm(length(mu))%*%R)
        if(fixed[k] || point.check(ts[k],x[1],x[2])) break
        ## If fail, return last fixed point
        if(r==100) return(k0)
      }
      xs[k,] <<- x

      ## Allow discontinuity at a fixed point
      if(fixed[k]) {
        x <- ms[k,] + drop(rnorm(4)%*%chol(Vs[,,k]))
        k0 <- k
      }
    }
    ## On success, return 0
    0
  }

  k <- n
  for(i in 1:50) {
    k <- if(i < 25) sample(k) else sample(n)
    if(k==0) return(data.frame(date=ts,lon=xs[,1],lat=xs[,2]))
  }
  NULL
}



##' Resample a Track by Linear Interpolation
##'
##' Linearly interpolate in longitude and latitude to resample a track
##' back to a regular time step.
##'
##' The input track is given as a dataframe where each row is an
##' observed location, with (at minimum) columns
##' \describe{
##' \item{'date'}{observation time (as GMT POSIXct),}
##' \item{'lon'}{observed longitude,}
##' \item{'lat'}{observed latitude}
##' }
##'
##' @title Interpolate Track
##' @param trk a dataframe representing a track
##' @param tstep the time step to reample to (in seconds)
##' @return a dataframe with columns
##'   \item{\code{date}}{observation time (as POSIXct)}
##'   \item{\code{lon}}{interpolated longitude}
##'   \item{\code{lat}}{interpolated latitude}
##' @importFrom stats approx
##' @export
interpolateTrack <- function(trk,tstep=60*60) {

  if(!is.null(trk)) {
    ts <- seq(min(trk$date),max(trk$date),tstep)
    data.frame(date=ts,
               lon=approx(as.numeric(trk$date),trk$lon,as.numeric(ts))$y,
               lat=approx(as.numeric(trk$date),trk$lat,as.numeric(ts))$y)

  }
}
