

makeTravisCoordsFromGRangesList <- function(comp, noBins) {
  
  # get all the check points
  tx_length <- as.numeric(sum(width(comp)))
  checkpoints_interval <- (tx_length-1)/(noBins-1)
  
  # get transcript name and strand
  tx_name <- names(comp)
  granges <- unlist(comp)
  tx <- granges[tx_name]
  strand <- as.character(as.data.frame(strand(tx))[[1]])
  chr <- as.character(as.data.frame(seqnames(tx))[[1]])
  
  # get coordinates
  t <- lapply(X=1:noBins, 
              FUN=.makeTravisCoordsForSingleIndex, 
              comp=comp,
              noBins=noBins,
              checkpoints_interval=checkpoints_interval,
              strand=strand,
              chr=chr,
              tx_name=tx_name)
 
  TravisCoords <- .combineListOfGRanges(t)
  
  # return the result
  return(TravisCoords)
}

.combineListOfGRanges <- function(t){
  txt <- "c(t[[1]]"
  for (i in 2:length(t)) {
    txt <- paste(txt,",t[[",i,"]]",sep="")
  }
  txt <- paste(txt,")",sep="")
  c <- parse(text=txt)
  
  # suppressWarnings
  # TravisCoords <- eval(c)
  TravisCoords <- suppressWarnings(eval(c))
  
  return(TravisCoords)
}
.makeTravisCoordsForSingleIndex <- function(
  index, comp, noBins, checkpoints_interval,
  strand, chr, tx_name) {
  
  # get checkpoints
  checkpoints <- round(1+checkpoints_interval*(index-1))
  checkpoints_transcript <- GRanges(seqnames=chr,
                                    IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                    strand=strand)
  
  # convert to genomic coordinates
  checkPoints_genomic <- pmapFromTranscripts(checkpoints_transcript, comp)
  
  # resize
  binWidth <- 4*round(checkpoints_interval)+1
  binWidth <- 31
  checkRegion_genomic <- resize(x=checkPoints_genomic, 
                                width=binWidth, 
                                fix="center")
  
  
  
  
  # add annotation information
  mcols(checkRegion_genomic) <- data.frame(txid=tx_name, pos=(index-1)/(noBins-1))
  
  return(checkRegion_genomic)
}