globalVariables(c(":="))

#Calculate Nominator 
#' Computation of the numerator of the moment estimator
#'
#' @param p1 numeric, frequencies in population 1
#' @param p2 numeric, frequencies in population 2
#' @param n1 numeric, number of gametes in population 1
#' @param n2 numeric, number of gametes in population 2
#'
#' @return a vector with the  numerators of the Fst moment estimator
#'
Compute_Nominator <- function (p1, p2, n1, n2){
  Nominator <- (p1-p2)^2-((p1*(1-p1))/(n1-1))-((p2*(1-p2))/(n2-1))  
}

#Calculate Nominator 
#' Computation of the numerator of the moment estimator
#'
#' @param p1 numeric, frequencies in population 1
#' @param p2 numeric, frequencies in population 2
#'
#' @return a vector with the denominators of the Fst moment estimator
#'
Compute_Denominator <- function(p1, p2){
  Denominator <- p1*(1-p2) + p2*(1-p1)
}

#Ratio of averages 
#' Computation of the Fst moment estimator
#'
#' @param Nominator numeric, numerators of the Fst moment estimator
#' @param Denominator numeric, denominators of the Fst moment estimator 
#'
#' @return a vector with the global Fst estimator
Ratio_Average <- function (Nominator,Denominator){
  if (is.matrix(Nominator)){
    Fst <- colMeans(Nominator)/colMeans(Denominator)
  } else {
    Fst <- mean(Nominator)/mean(Denominator)
  }
  return(Fst)
}



#' Convert the Freq and NbGamete tables into a list.  
#'
#' @param Freq A data.frame or matrix of frequencies where each row corresponds to a marker, each column corresponds to a population,
#' @param NbGamete A data.frame or matrix of number of gametes where each row corresponds to a marker, and each column corresponds to a population
#'
#' @return a list of data.frames, each corresponding to a population. 
#' @export
#'
#' @description The function builds a list where each element corresponds to a population present in both Freq and NbGametes (all other
#' populations being discarded). Each element consists of a data.frame with 2 columns, Freq and NbGamete.
#' @examples 
#' ## Load the HGDP data 
#' data(Freq);data(NbGamete)
#' FreqNbG <- BuildFreqNbG(Freq,NbGamete)
BuildFreqNbG <- function(Freq,NbGamete){
  
  Freq <- as.data.frame(Freq)
  NbGamete <- as.data.frame(NbGamete)
  
  ## Check names 
  if(!identical(names(Freq),names(NbGamete))){
    message("Populations in Freq and NbGamete are different...")
  }
  Population.list <- intersect(names(Freq),names(NbGamete)) %>% setNames(.,.)
  
  ## Check Nb Markers
  if(nrow(Freq)!=nrow(NbGamete)){
    stop("The numbers of rows (markers) in Freq and NbGamete are different")
  }
  
  ## Build the FreqNbG object
  FreqNbG <- map(Population.list, ~data.frame(Freq=Freq[[.x]],NbGamete=NbGamete[[.x]]))
  return(FreqNbG)
}


#' HudsonFst.m
#'
#' Compute Hudson Fst moment estimator at marker level
#'
#' @param FreqNbG a list of data.frames (one per population) with 2 columns: Freq and NbGamete  
#'
#' @return a list of data.frames with 2 columns: Fst and Weight.
#' @export
#' @examples 
#' ## Load the FreqNbG object build from the HGDP data
#' data(Freq);data(NbGamete)
#' FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' HFst.m <- HudsonFst.m(FreqNbG)
HudsonFst.m <- function(FreqNbG){
  
  ## Basic quantities
  NbSnp <- nrow(FreqNbG[[1]])
  NbPop <- length(FreqNbG)

  ## Get the pairs for which Fst is required
  PopNames <- names(FreqNbG)
  PairPop <- combn(1:NbPop , 2 ,simplify = F)
  PairPopNames <- sapply(PairPop, function(pair){
    paste0(PopNames[pair[1]],'_',PopNames[pair[2]])
  })
  
  #Compute num and den
  Res <- lapply(PairPop, FUN= function(pair){
    Num <- Compute_Nominator(FreqNbG[[pair[1]]]$Freq, FreqNbG[[pair[2]]]$Freq, 
                             FreqNbG[[pair[1]]]$NbGamete, FreqNbG[[pair[2]]]$NbGamete)
    Den <- Compute_Denominator(FreqNbG[[pair[1]]]$Freq, FreqNbG[[pair[2]]]$Freq)
    Fst <- Num/Den
    Fst[which(is.na(Fst))] <- 0
    return(data.frame(Fst=Fst,Weight=Den))
  })
  names(Res) <- PairPopNames
  return(Res)  
}





#' HudsonFst.gw
#' 
#' Compute genome-wide Hudson Fst moment estimator
#'
#' @param HFst.m a list of data.frames as obtained with function \code{HudsonFst.m}
#' @param Mat boolean, should the result be output as a matrix. 
#'
#' @return By default a matrix of Hudson Fst coefficients, a vector otherwise.
#' @import stringr
#' @export
#' @examples 
#' data(Freq);data(NbGamete)
#' FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' HFst.m <- HudsonFst.m(FreqNbG)
#' HFst.chr <- HudsonFst.gw(HFst.m)
HudsonFst.gw <- function(HFst.m,Mat=TRUE){
  
  ## Compute the genome (chromosome in practice) wide Fst for each pair
  Res <- sapply(HFst.m, function(dt){
    sum(dt$Fst*dt$Weight)/sum(dt$Weight)
  })
  
  ## Display as a matrix if required
  if (Mat){
    ## Get the pop names
    PairNames <- names(HFst.m) %>% 
      str_split(pattern='_')
    PopNames <- PairNames %>% 
      map(~.[2]) %>% 
      unique %>% 
      as.character %>% 
      c(PairNames[[1]][1],.)
    NbPops <- length(PopNames)
    ## Build Fst matrix
    FstMatrix <- matrix(0,NbPops,NbPops)
    FstMatrix[lower.tri(FstMatrix)] <- Res
    FstMatrix[upper.tri(FstMatrix)] <- t(FstMatrix)[upper.tri(FstMatrix)]
    colnames(FstMatrix) <- row.names(FstMatrix) <- PopNames
    Res <- FstMatrix
  }
  return(Res)
}


#' Summarise1Profile
#'
#' @param profile a vector, corresponding to the Fst profile of a pair of populations
#' @param snpinfo a data.frame, providing information about markers
#' @import dplyr tibble purrr
#' @importFrom rlang .data
#' @return a data.frame that combines both objects
Summarise1Profile <- function(profile,snpinfo){
  Segment <- profile %>% 
    diff %>% 
    `!=`(0) %>% 
    which %>% c(0,.,length(profile)) %>% 
    diff %>% 
    rep(1:length(.),times=.)  
  DF <- data.frame(snpinfo,profile,Segment) %>%
    as_tibble %>% 
    group_by(Segment) %>% 
    summarise(Start=min(.data$POS),End=max(.data$POS),NbSnp=n(),Fst=mean(profile),.groups='drop') 
  return(DF)
}

#' ProfilingSummary
#' 
#' Summary of Fst profiles
#'
#' @param FstProfiles a list of Fst profiles as obtained from function \code{HudsonFst.prof} 
#' @param SnpInfo a data.frame providing information about markers
#'
#' @return a list of data.frame. Each data.frame summarizes a Fst profile, in terms of number of segments, 
#' Start and End positions, length (i.e. number of markers)  and Fst level of each segment.
#' @export
#' 
#' @examples 
#' data(Freq);data(NbGamete);data(Info)
#' FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' HFst.m <- HudsonFst.m(FreqNbG)
#' TwoPops <- list(First="Colombian",Second="Tuscan")
#' HFst.prof <- HudsonFst.prof(HFst.m,Contrast=TwoPops)
#' PS <- ProfilingSummary(HFst.prof,Info)
ProfilingSummary <- function(FstProfiles,SnpInfo){
  Res <- map(FstProfiles,Summarise1Profile,snpinfo=SnpInfo) 
  return(Res)
}


#' Filtering markers based on allelic frequencies
#'
#' @param FreqNbG a list of data.frames (one per population) with 2 columns: Freq and NbGamete
#' @param Maf a numerci value for the thresholding of minor allelic frequencies 
#'
#' @return a vector of positions to be removed
#' @export
Freq.filt <- function(FreqNbG,Maf=0){

  Nb1 <- map(FreqNbG, ~ .x$Freq*.x$NbGamete) %>% 
    bind_cols %>% 
    rowSums
  N <- map(FreqNbG, ~ .x$NbGamete) %>% 
    bind_cols %>% 
    rowSums
  Rmv.idx <- (((Nb1/N)>= 1-Maf)|((Nb1/N)<= Maf)) %>% 
    which
  return(Rmv.idx)
}


#' @title HudsonFst.prof
#'
#' @description Perform FST profiling between pairs of pops, as requested by Contrast. If no contrast is provided, all pairs are considered
#'
#' @param HFst.m A list of data.frame with two columns each, Fst and Weight, 
#' as provided by the \code{HudsonFst.m} function
#' @param Contrast a list of two vectors with the names of the populations to be contrasted 
#' @param Kmax maximum number of breakpoints to be considered
#' @param NbSegCrit the criterion used for the choice of the number of segments 
#' @param parallel a boolean, should the profiling be parallelized (using future) or not
#' @return a smoothed profile
#' @import furrr future
#' @export
#' 
#' @examples 
#' data(Freq);data(NbGamete)
#' FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' HFst.m <- HudsonFst.m(FreqNbG)
#' 
#' ## Two population analysis
#' TwoPops <- list(First="Colombian",Second="Tuscan")
#' HFst.prof <- HudsonFst.prof(HFst.m,Contrast=TwoPops)
#' 
#' ## The full example execution takes a few seconds.
#' ## Two sets of populations to contrast
#' \donttest{
#' Contrast <- list(America=c("Colombian","Maya"),Europe=c("Tuscan","Italian"))
#' Profiles <- HudsonFst.prof(HFst.m,Contrast=Contrast)
#' }
#' ## For larger lists and/or larger marker sets, 
#' ## use the future package for parallel computation:
#' \donttest{
#' future::plan("multisession",workers=4)
#' Profiles <- HudsonFst.prof(HFst.m,Contrast=Contrast)
#' future::plan("default")
#' }
HudsonFst.prof <- function(HFst.m,Contrast=NULL,Kmax=100,NbSegCrit='biggest.S3IB',parallel=TRUE){
  
  ## Checkings
  if(!is.list(HFst.m)){
    stop('HFst.m should be a list')
  }
  if(!(is.list(Contrast)|is.null(Contrast))){
    stop('Contrast should be a list')
  }
  
  ## Get the list of pairwise profiles to be computed
  HFstNames <- names(HFst.m)
  if(is.null(Contrast)){
    PairNames <- HFstNames
  } else if (length(Contrast)==1){
    InvolvedPairs <- combn(Contrast[[1]],2) %>% t 
    InvolvedPairs <- c(apply(InvolvedPairs,1,paste,collapse='_'),apply(InvolvedPairs[,2:1],1,paste,collapse='_'))
    PairNames <- map_lgl(HFstNames,~ .x %in% InvolvedPairs) %>% which %>% HFstNames[.]
  } else {
    InvolvedPairs <- c(outer(Contrast[[1]],Contrast[[2]],paste,sep='_') %>% as.character,
                       outer(Contrast[[2]],Contrast[[1]],paste,sep='_') %>% as.character)
    PairNames <- map_lgl(HFstNames,~ .x %in% InvolvedPairs) %>% which %>% HFstNames[.]
  }
  
  ## Perform Fst profiling
  if(parallel){
    Profiles <- furrr::future_map(HFst.m[PairNames],MakeProfile.sn_nomemory,Kmax=Kmax,method=NbSegCrit)
  } else {
    Profiles <- map(HFst.m[PairNames],MakeProfile.sn_nomemory,Kmax=Kmax,method=NbSegCrit)
  }
  return(Profiles)    
}





#' DF4Plot1Prof
#' 
#' Shape the data for Fst profile representation
#'
#' @param Info a data.frame providing information about markers
#' @param HF a data.frame with 2 columns Fst and Weight, as obtained from the \code{HudsonFst.m} function  
#' @param FstProf an Fst profile, as obtained from the \code{HudsonFst.prof} function  
#' @param Coord a vector with the minimum and maximum coordinates (i.e. positions along the genome), 
#' providing the range of the genomic region that will be plotted, optional. 
#' @param Threshold a numeric value. Markers belonging to regions whose Fst profile is higher than threshold
#' will be highlighted. Optional.
#'
#' @importFrom rlang .data
#' @return a data.frame that can be used as an input for function \code{Plot1Prof}
DF4Plot1Prof <- function(Info,HF,FstProf,Coord=NULL,Threshold=NULL){
  if(is.null(Threshold)){
    Threshold=10 # A value that cannot be reached !
  }
  if(is.null(Coord)){
    Coord=range(Info$POS)
  }
  DF <- data.frame(Info,HF,FstProf) %>% 
    mutate(Detected = (FstProf>=Threshold)) %>% 
    filter(.data$POS>=Coord[1] & .data$POS<=Coord[2]) %>% 
    dplyr::select(contains('CHR'),contains('Chr'),.data$POS,.data$Weight,.data$Fst,FstProf,.data$Detected) %>% 
    as_tibble
  return(DF)
}



#' ContrastSummary
#' 
#' Summarize multiple Fst profiles 
#'
#' @param PS a list of profile summaries, as provided by the \code{ProfilingSummary} function
#' @param RefLevel a list of reference (i.e. baseline) Fst levels
#' @param Ratio.thres a numeric value, regions exhibiting Fst levels whose ratio with the reference level is 
#' higher than Ratio.thres will be highlighted.
#' @param NbSnp.min an integer. The minimum number of markers required to highlight a region 
#'
#' @return a tibble
#' @importFrom rlang .data 
#' @export
#' @examples 
#' ## The full example execution takes a few seconds.
#' # data(Freq);data(NbGamete)
#' # FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' # HFst.m <- HudsonFst.m(FreqNbG)
#' 
#' ## Two sets of populations to contrast
#' # Contrast <- list(America=c("Colombian","Maya"),Europe=c("Tuscan","Italian"))
#' # Profiles <- HudsonFst.prof(HFst.m,Contrast=Contrast)
#' # PS <- ProfilingSummary(Profiles,Info) 
#' 
#' # RefLevel <- rapply(Profiles,median,classes = "numeric",how='list')
#' # Ratio.thres <- 3
#' # NbSnp.min <- 1
#' # CS <- ContrastSummary(PS, RefLevel,
#' #                       Ratio.thres=Ratio.thres,
#' #                       NbSnp.min=NbSnp.min)
ContrastSummary <- function(PS, RefLevel, Ratio.thres=3,NbSnp.min=1){

  Res <- map2(PS,RefLevel, ~{
    mutate(.x,RefLevel=.y,
           Ratio=Fst/RefLevel,
           Select = (Ratio >Ratio.thres & NbSnp > NbSnp.min)) %>% 
      {tibble(NbSel=rep(.$Select, .$NbSnp) %>% as.numeric,
              MeanRatio=rep(.$Ratio, .$NbSnp))}
    }) %>% 
    Reduce('+',.) %>% 
    mutate(MeanRatio=.data$MeanRatio/length(PS)) %>% 
    as_tibble
  return(Res)
}


#' MergeRegion
#' 
#' Merge adjacent top regions
#'
#' @param DT a data.frame 
#' @param Crit a string corresponding to the name of the criterion used for selecting the top regions
#'
#' @return a simplified data.frame (tibble)
MergeRegion <- function(DT,Crit){
  Try=TRUE
  while(Try){
    Need <- diff(DT$Segment)  
    if(any(Need==1)){
      Who <- which(Need==1)[1]
      DT$Start[Who+1] <- DT$Start[Who]
      DT$NbSnp[Who+1] <- DT$NbSnp[Who]+DT$NbSnp[Who+1]
      DT[[Crit]][Who+1] <- min(DT[[Crit]][Who],DT[[Crit]][Who+1])
      DT <- DT[-Who,]
    } else {
      Try=FALSE  
    }
  }
  return(DT)
}


#' ContrastTopRegions
#'
#' @param CS a list of contrast summaries as obtained from function \code{ContrastSummary}
#' @param Crit a string providing the name of the variable to use to select regions
#' @param Info a data.frame providing information about markers
#' @param Thres the threshold to be used on the Crit variable
#' @param Simplify a boolean specifying whether the results should be displayed 
#' as a list (by-default option) or as a single data.frame
#'
#' @return a data.frame or a list of data.frames
#' @importFrom rlang .data
#' @import tidyr 
#' @export
#' @examples 
#' ## The full example execution takes a few seconds.
#' # data(Freq);data(NbGamete)
#' # FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' # HFst.m <- HudsonFst.m(FreqNbG)
#' 
#' ## Two sets of populations to contrast
#' # Contrast <- list(America=c("Colombian","Maya"),Europe=c("Tuscan","Italian"))
#' # Profiles <- HudsonFst.prof(HFst.m,Contrast=Contrast)
#' # PS <- ProfilingSummary(Profiles,Info) 
#' 
#' # RefLevel <- rapply(Profiles,median,classes = "numeric",how='list')
#' # Ratio.thres <- 3
#' # NbSnp.min <- 1
#' # CS <- ContrastSummary(PS, RefLevel,
#' #                       Ratio.thres=Ratio.thres,
#' #                       NbSnp.min=NbSnp.min)
#' # NbSel.thres <- 2
#' # TopRegions <- ContrastTopRegions(CS = CS,Crit = 'NbSel',Info = Info,
#' #                                 Thres = NbSel.thres, Simplify=TRUE)
ContrastTopRegions <- function(CS,Crit,Info,Thres,Simplify=FALSE){
  if(is.data.frame(CS)){
    if(!is.data.frame(Info)){
      stop("CS is a data.frame, so Info should be a data.frame too.")
    }
    CS <- list(CS)
    Info <- list(Info)
  }
  TopRegions <- map(CS,~.x[[Crit]]) %>% 
    map2(.,Info,Summarise1Profile) %>% 
    imap(~mutate(.x,Chromosome=.y)) %>% 
    bind_rows %>%
    rename(!!Crit:=.data$Fst) %>% 
    {.[.[[Crit]]>=Thres,]} %>% 
    arrange(.data$Chromosome,.data$Segment)
  if(Simplify){
    TopRegions <- TopRegions %>% 
      tidyr::nest(data=-.data$Chromosome) %>% 
      mutate(data=map(data,MergeRegion,Crit=Crit)) %>% 
      unnest(cols=c(data)) %>% 
      arrange(.data$Chromosome)
  }
  return(TopRegions)
}

