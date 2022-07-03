
#' RawPlot
#' 
#' Plot the Fst estimates along a (portion of) chromosome
#'
#' @param Info a data.frame providing information about markers
#' @param HF a data.frame with 2 columns, Fst and Weight, as provided by the \code{HudsonFst.m} function
#' @param Coord a vector with the minimum and maximum coordinates (i.e. positions along the genome), 
#' providing the range of the genomic region that will be plotted. 
#' @param Title a string providing a title for the graph.
#' @import ggplot2  
#' @importFrom scales unit_format
#' @importFrom rlang .data
#' @return a ggplot object.
#' @export
RawPlot <- function(Info,HF,Coord=NULL,Title=''){
  if(is.null(Coord)){
    Coord=range(Info$POS)
  }
  DF <- data.frame(Info,HF) %>% 
    filter(.data$POS>=Coord[1] & .data$POS<=Coord[2]) %>% 
    dplyr::select(contains('CHR'),contains('Chr'),.data$POS,.data$Weight,.data$Fst) %>%
    mutate(Weight_10=.data$Weight/10) %>% 
    as_tibble
  Range = range(DF$Fst)
  if((Coord[2]-Coord[1])/(1e6)>2){
    Unit = "Mb"; Scale = 1e-6
  } else {
    Unit = "Kb"; Scale = 1e-3
  }
  Labs <- labs(title = Title,x=NULL,y=NULL)
  GG <- ggplot() +
    scale_x_continuous(labels = scales::unit_format(unit = Unit, scale = Scale)) +
    scale_size(range=c(0.1,2),guide = "none") +
    scale_color_manual(breaks = c("0", "1"), values=c("black","blue"),guide="none") +
    theme(axis.text.x = element_text(angle=45,size=12),
          axis.title=element_text(size=15),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=15),
          plot.margin = margin(1,2,-2,1,"mm")) +
    ylim(Range) +
    geom_point(data=DF, aes_string(x = 'POS',y='Fst',size='Weight_10')) +
    Labs
  return(GG)
}


#' Plot1Prof
#' 
#' Display the graphical representation of an Fst profile
#'
#' @param DF a data.frame, as provided by function \code{DF4Plot1Prof}
#' @param Title a string providing a title for the graph, optional. 
#' @param Range a vector with the minimum and maximum values for the y-axis 
#'
#' @return a ggplot object
#' @importFrom rlang .data
Plot1Prof <- function(DF,Title='',Range=NULL){
  DF <- DF %>% 
    mutate(PointColor = as.character(as.numeric(.data$Detected)), LineSize=0.5) %>% 
    mutate(Weight_10=.data$Weight/10)
  if (is.null(Range)){
    Range <- range(DF$Fst)
  }
  PosRange <- range(DF$POS)
  if((PosRange[2]-PosRange[1])/(1e6)>2){
    Unit = "Mb"; Scale = 1e-6
  } else {
    Unit = "Kb"; Scale = 1e-3
  }
  Labs <- labs(title = Title,x=NULL,y=NULL)
  GG <- ggplot() +
    scale_x_continuous(labels = unit_format(unit = Unit, scale = Scale)) +
    scale_size(range=c(0.1,2),guide = "none") +
    scale_color_manual(breaks = c("0", "1"), values=c("black","blue"),guide="none") +
    theme(axis.text.x = element_text(angle=45,size=12),
          axis.title=element_text(size=15),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=15),
          plot.margin = margin(1,2,-2,1,"mm")) +
    ylim(Range) +
    geom_point(data=DF, aes_string(x = 'POS',y='Fst',colour='PointColor',size='Weight_10')) +
    geom_line(aes_string(x='POS',y='FstProf'),color='red',size=1.5,data=DF) +
    Labs
  return(GG)
}



#' ContrastGraphSummary
#' 
#' Display mean ratio and/or number of selection graphs
#'
#' @param CS a contrast summary, as provided by function \code{ContrastSummary}
#' @param Info a data.frame providing information about markers
#' @param Ratio.thres a numeric value, regions exhibiting Fst levels whose ratio with the reference level is 
#' higher than Ratio.thres will be highlighted. 
#' @param Coef a scalar, controling font sizes for the graph, optional 
#' @param CutNbSel a scalar, providing a y-value for an horizontal line on the NbSel graph 
#' @param CutMeanRatio a scalar, providing a y-value for an horizontal line on the MeanRatio graph
#' @importFrom stats lm mad pt setNames weighted.mean
#' @return two ggplots objects, called NbSel and MeanRatio, respectively.
#' @export
ContrastGraphSummary <- function (CS, Info, Ratio.thres, Coef = 1, CutNbSel = NULL, CutMeanRatio = NULL){
  TH <- theme(axis.text.x = element_text(angle = 45, size = 12), 
              axis.title = element_text(size = 15), axis.text.y = element_text(size = 10), 
              plot.title = element_text(size = 15), plot.margin = margin(1, 
                                                                         2, -2, 1, "mm"))
  Tmp <- bind_cols(Info, CS)
  PosRange <- range(Tmp$POS)
  if ((PosRange[2] - PosRange[1])/(1e+06) > 2) {
    Unit = "Mb"
    Scale = 1e-06
  } else {
    Unit = "Kb"
    Scale = 0.001
  }
  GG.NbSel <- ggplot(Tmp, aes_string(x = 'POS', y = 'NbSel')) + 
    scale_x_continuous(labels = unit_format(unit = Unit, scale = Scale)) + 
    scale_size(range = c(0.1, 2), guide = "none") + 
    scale_color_manual(breaks = c("0", "1"), 
                       values = c("black", "blue"), guide = "none") + 
    geom_line() + ylab(paste0("Nb prof. > ", Ratio.thres)) + 
    xlab("") + TH + ggtitle("")
  if (!is.null(CutNbSel)) {
    GG.NbSel <- GG.NbSel + geom_hline(yintercept = CutNbSel, 
                                      color = "red", size = 1.2)
  }
  GG.MeanRatio <- ggplot(Tmp, aes_string(x = 'POS', y = 'MeanRatio')) + 
    scale_x_continuous(labels = unit_format(unit = Unit, scale = Scale)) + 
    scale_size(range = c(0.1, 2), guide = "none") + 
    scale_color_manual(breaks = c("0", "1"),
                       values = c("black", "blue"), guide = "none") + 
    geom_line() + ylab("Mean Ratio") + 
    xlab("") + TH + ggtitle("")
  if (!is.null(CutMeanRatio)) {
    GG.MeanRatio <- GG.MeanRatio + geom_hline(yintercept = CutMeanRatio, 
                                              color = "red", size = 1.2)
  }
  return(list(NbSel = GG.NbSel, MeanRatio = GG.MeanRatio))
  
}




#' HeatMap
#' 
#' Make a frequency heatmap 
#'
#' @param Min the starting position value of the region
#' @param Max the end position value of the region
#' @param chr a string providing the chromosome name, optional
#' @param Info a data.frame providing information about markers
#' @param FreqNbG a list of data.frames (one per population) with two columns: Freq and NbGamete
#' @param Dir a string providing the name of the directory where the graph should be saved, optional
#' @param Weights a vector of weights associated with each marker, optional
#' @param Weight.thres a numeric value. Markers with weights lower than this threshold
#' will be discarded from the graphical representation. Optional.
#' @param NbAdjM an integer providing the number of markers before and after the highlighted regions 
#' that should be added to the graphical representation, optional.  
#' @param Subsets a list of character vectors with the population names, optional.
#' @return A heatmap where rows correspond to markers, and columns to populations. 
#' @import gplots ggplot2
#' @importFrom grDevices dev.off png
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
#' #                                  Thres = NbSel.thres, Simplify=TRUE)
#' # HeatMap(Min = TopRegions[1,]$Start,
#' #         Max = TopRegions[1,]$End,
#' #         chr = TopRegions[1,]$Chromosome,
#' #         Info = Info,
#' #         FreqNbG = FreqNbG,
#' #         Subsets = Contrast)                                
HeatMap <- function(Min,Max,chr=NULL,Info,FreqNbG,Dir=NULL,Weights=NULL,Weight.thres=0.05,NbAdjM=0,Subsets=NULL){
  
  ## Get the markers of interest
  Cand <- which((Info$POS>=Min) & (Info$POS<=Max))
  if(!is.null(Weights)){
    Cand <- Weights[Cand] %>% 
      `>`(Weight.thres) %>% 
      which %>% 
      Cand[.]
  }
  
  ## Add some adjacent markers if required
  if(NbAdjM>0){
    RSC <- rep(c("white", "black","white"), times =
                 c(NbAdjM,length(Cand),NbAdjM))
    Cand <- c((min(Cand)-NbAdjM):(min(Cand)-1),
              Cand,
              (max(Cand)+1):(max(Cand)+NbAdjM))
  } else {
    RSC <- rep("white",length(Cand))
  }
  
  ## Reshape the frequency dataset
  freq <- FreqNbG %>% map(~.x$Freq) %>% bind_cols() %>% as.matrix
  #If subsets are specified
  if(is.null(Subsets)){
    ColSep=ncol(freq)
  } else {
    Subsets.all <- Subsets %>% Reduce('c',.)
    NbPerSub <- map_int(Subsets,length)
    #Check that all pops are present
    if (any(!(Subsets.all %in% colnames(freq)))) {
      stop('Some pops in Subsets are not described in Freq')
    }
    freq <- map(Subsets, ~ freq %>% as.data.frame %>% select(one_of(.x))) %>% bind_cols()
    ColSep=cumsum(NbPerSub)
  }
  row.names(freq) <- Info$MARKER
  
  ## If require build the .png file name
  if(!is.null(Dir)){
    if(is.null(chr)){
      FN <- paste0('Heatmap_',Min,'to',Max)  
    } else {
      FN <- paste0('HeatmapChr',chr,'_',Min,'to',Max)
    }
    if(!is.null(Weights)){
      FN <- paste0(FN,'_weighted')
    }
    png(paste0(Dir,FN,'.png'))
  }
  
  ## Make the heatmap
  
  freq[Cand,] %>% 
    as.matrix %>% 
    heatmap.2(., col=greenred(75),
              breaks=seq(0,1,length.out=76),
              Rowv = NULL,Colv = NULL,
              trace = 'none',
              dendrogram = 'none',
              colsep=c(ColSep),
              # keysize = 0.9,
              key.par = list(cex=0.8),
              key.title = '',
              key.xlab = '',
              key.ylab = '',
              density.info="none",
              cexCol =  2,
              margins=c(9.5,6),
              RowSideColors = RSC)
  
  ## Close file (if needed)
  if(!is.null(Dir)){
    dev.off()
  }
}

# Dir=NULL
# Weights=NULL
# Weight.thres=0.05
# NbAdjM=0
# Min = TopRegions[1]$Start
# Max = TopRegions[1]$End
# chr = TopRegions[1]$Chromosome
# Info = Info
# FreqNbG = FreqNbG
# Subsets = Contrast


#' Plot Fst values along chromosomes
#'
#' @param Info a data.frame providing information about markers
#' @param HFst.m a data.frame with 2 columns, Fst and Weight, as provided by the \code{HudsonFst.m} function  
#' @param HFst.prof a data.frame corresponding to one item of the output of the \code{HudsonFst.prof} function
#' @param Coord a vector with the minimum and maximum coordinates (i.e. positions along the genome)
#' providing the range of the genomic region that will be plotted. 
#' @param Ref a value to plot a reference line 
#' @param Threshold a value to plot a threshold line
#'
#' @return a ggplot object
#' @import ggplot2  
#' @export
#'
#' @examples 
#' ## The full example execution takes a few seconds.
#' data(Freq);data(NbGamete)
#' FreqNbG <- BuildFreqNbG(Freq,NbGamete)
#' HFst.m <- HudsonFst.m(FreqNbG)
#' TwoPops <- list(First="Colombian",Second="Tuscan")
#' HFst.prof <- HudsonFst.prof(HFst.m,Contrast=TwoPops)
#' 
#' ## Plot the raw data
#' \donttest{
#' HudsonFst.plot(Info,HFst.m$Colombian_Tuscan)
#' }
#' ## Plot the raw data and the segmentation
#' \donttest{
#' HudsonFst.plot(Info,HFst.m$Colombian_Tuscan,HFst.prof$Colombian_Tuscan)
#' }
#' ## Add a background/reference level
#' \donttest{
#' RefLevel <- median(HFst.prof$Colombian_Tuscan)
#' HudsonFst.plot(Info,HFst.m$Colombian_Tuscan,HFst.prof$Colombian_Tuscan,
#'                Ref=RefLevel)
#' }
#' ## Add a threshold 
#' \donttest{
#' Threshold <- 3*RefLevel
#' HudsonFst.plot(Info,HFst.m$Colombian_Tuscan,HFst.prof$Colombian_Tuscan,
#'                Ref=RefLevel,Threshold = Threshold)
#' }
HudsonFst.plot <- function(Info,HFst.m,HFst.prof=NULL,Coord=NULL,Ref=NULL,Threshold=NULL){
  
  if(is.null(HFst.prof)){
    ## Plot the per marker Fst estimates
    GG <- RawPlot(Info,HFst.m)
  } else{
    ## Plot marker estimates + profiles
    GG <- DF4Plot1Prof(Info = Info, HF = HFst.m, FstProf = HFst.prof) %>% 
      Plot1Prof
  }
  if (!is.null(Ref)){
    ## Add a reference line
    GG <- GG + geom_hline(yintercept = Ref,col='green',size=1.5)
  }
  if (!is.null(Ref)){
    ## Add a threshold line
    GG <- GG + geom_hline(yintercept = Threshold,col='blue',size=1.5)
  }
  
  return(GG)
}

