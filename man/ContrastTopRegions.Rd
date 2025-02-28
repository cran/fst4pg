% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FstBattia.R
\name{ContrastTopRegions}
\alias{ContrastTopRegions}
\title{ContrastTopRegions}
\usage{
ContrastTopRegions(CS, Crit, Info, Thres, Simplify = FALSE)
}
\arguments{
\item{CS}{a list of contrast summaries as obtained from function \code{ContrastSummary}}

\item{Crit}{a string providing the name of the variable to use to select regions}

\item{Info}{a data.frame providing information about markers}

\item{Thres}{the threshold to be used on the Crit variable}

\item{Simplify}{a boolean specifying whether the results should be displayed 
as a list (by-default option) or as a single data.frame}
}
\value{
a data.frame or a list of data.frames
}
\description{
ContrastTopRegions
}
\examples{
## The full example execution takes a few seconds.
# data(Freq);data(NbGamete)
# FreqNbG <- BuildFreqNbG(Freq,NbGamete)
# HFst.m <- HudsonFst.m(FreqNbG)

## Two sets of populations to contrast
# Contrast <- list(America=c("Colombian","Maya"),Europe=c("Tuscan","Italian"))
# Profiles <- HudsonFst.prof(HFst.m,Contrast=Contrast)
# PS <- ProfilingSummary(Profiles,Info) 

# RefLevel <- rapply(Profiles,median,classes = "numeric",how='list')
# Ratio.thres <- 3
# NbSnp.min <- 1
# CS <- ContrastSummary(PS, RefLevel,
#                       Ratio.thres=Ratio.thres,
#                       NbSnp.min=NbSnp.min)
# NbSel.thres <- 2
# TopRegions <- ContrastTopRegions(CS = CS,Crit = 'NbSel',Info = Info,
#                                 Thres = NbSel.thres, Simplify=TRUE)
}
