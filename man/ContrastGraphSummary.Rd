% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Graphs.R
\name{ContrastGraphSummary}
\alias{ContrastGraphSummary}
\title{ContrastGraphSummary}
\usage{
ContrastGraphSummary(
  CS,
  Info,
  Ratio.thres,
  Coef = 1,
  CutNbSel = NULL,
  CutMeanRatio = NULL
)
}
\arguments{
\item{CS}{a contrast summary, as provided by function \code{ContrastSummary}}

\item{Info}{a data.frame providing information about markers}

\item{Ratio.thres}{a numeric value, regions exhibiting Fst levels whose ratio with the reference level is 
higher than Ratio.thres will be highlighted.}

\item{Coef}{a scalar, controling font sizes for the graph, optional}

\item{CutNbSel}{a scalar, providing a y-value for an horizontal line on the NbSel graph}

\item{CutMeanRatio}{a scalar, providing a y-value for an horizontal line on the MeanRatio graph}
}
\value{
two ggplots objects, called NbSel and MeanRatio, respectively.
}
\description{
Display mean ratio and/or number of selection graphs
}
