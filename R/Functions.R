#' @import utils
## quiets concerns of R CMD check re: the .'s that appear in pipelines
#if(getRversion() >= "2.15.1")  utils::globalVariables(c(".",">"))
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))




################################################################################

#' @title MakeProfile.op
#'
#' @description Perform segmentation on a given dataset and returns the segmented profile
#' @param DF a data.frame with two columns, Fst and Weights, as provided by the \code{HudsonFst.m} function
#' @param coef.pen.value coef to use for penaly 2*coef.pen.value*log(n)
#' @param sd.y a numeric value corresponding to the (estimated) standard deviation of the signal. If NULL (default) the value is automatically estimated.
#' @return a smoothed profile
#' @importFrom graphics abline legend lines par plot
#' @export
MakeProfile.op <- function(DF, coef.pen.value=1,sd.y=NULL){

  ## Checks
  if(!(is.null(sd.y)|is.numeric(sd.y))){
    stop("sd.y should be a numeric value")
  }
  Y <- DF$Fst
  Weights <- DF$Weight
  
  IdxNotZero <- which(Weights > 0)
  Y <- Y[IdxNotZero]
  Weights <- (Weights[IdxNotZero])
  n <- length(Y)
  
  ## Standardize signal
  if(is.null(sd.y)){
    sd.y <- fpopw::sdDiff(Y)
  } 
  Y <- (Y/sd.y)
  
  ## recover only contrast with Fpop_w
  res.op <- fpopw::Fpop_w(Y, Weights, lambda=2*coef.pen.value*log(n))
  tauHat <- res.op$t.est
  tauHat.comp <- IdxNotZero[tauHat]
  if(tauHat.comp[length(tauHat.comp)]!=length(DF$Fst)){
      tauHat.comp[length(tauHat.comp)] <- length(DF$Fst)
    }

  smt <- fpopw::getSMT_(DF$Fst, DF$Weight, tauHat.comp)

}

################################################################################

#' @title MakeProfile.sn_nomemory
#'
#' @description Perform segmentation on a given dataset and returns the segmented profile
#' @param DF a data.frame with two columns, Fst and Weights, as provided by the \code{HudsonFst.m} function
#' @param Kmax max number of changes used for model selection (check with crops)
#' @param method a string, the name of the criterion used for model selection: (1) "givenVariance" = using the penalty of Lebarbier 2005 given a estimator of the variance, (2) "biggest.S3IB" = biggest=TRUE in saut taken from S3IB, (3) "notbiggest.S3IB"
#' To be chosen amongs3ib.jump, s3ib.nojump, pre, ddse (capushe) or jump (capushe)
#' @param sd.y a numeric value corresponding to the (estimated) standard deviation of the signal. If NULL (default) the value is automatically estimated.
#' @return a smoothed profile
#' @importFrom graphics abline legend lines par plot
#' @export
MakeProfile.sn_nomemory <- function(DF, Kmax , method='biggest.S3IB',sd.y=NULL){

  ## Checks
  if(!(is.null(sd.y)|is.numeric(sd.y))){
    stop("sd.y should be a numeric value")
  }
  
  Y <- DF$Fst
  Weights <- DF$Weight
  
  IdxNotZero <- which(Weights > 0)
  Y <- Y[IdxNotZero]
  Weights <- (Weights[IdxNotZero])
  n <- length(Y)
  
  ## Standardize signal
  ## Standardize signal
  if(is.null(sd.y)){
    sd.y <- fpopw::sdDiff(Y)
  } 
  Y <- (Y/sd.y)

  ## recover only contrast with Fpsn no memory
  res.sn <- fpopw::Fpsn_w_nomemory(Y, Weights, Kmax=Kmax)
  K_sel <- fpopw::select_Fpsn(res.sn, method=method,sigma=1) 
  tauHat <- fpopw::getTau_nomemory(res.sn, K_sel)
  ## this should be replace using simpleSelection from fpopw...{
  ## select K
  tauHat.comp <- IdxNotZero[tauHat]
#  if(tauHat.comp[length(tauHat.comp)] != length(DF$Fst)) tauHat.comp <- c(tauHat.comp, length(DF$Fst))
  if(tauHat.comp[length(tauHat.comp)]!=length(DF$Fst)){
      tauHat.comp[length(tauHat.comp)] <- length(DF$Fst)
    }
  
  smt <- fpopw::getSMT_(DF$Fst, DF$Weight, tauHat.comp)

}


 


