

TravisPlot <- function(peak,TravisCoordsFromTxDb=NA,txdb=NA,genome=NA,saveToPDFprefix=NA){
  
  # make sure the travis coordinates are available
  if (is.na(TravisCoordsFromTxDb)&is.na(txdb)&is.na(genome)) {
    stop("Most provide one of the three: TravisCoords, txdb or genome")
  } 
  
  if ( suppressWarnings(is.na(TravisCoordsFromTxDb)) ) {
    if (suppressWarnings(is.na(txdb))) {
      print("Downloading Transcriptome Information from UCSC ...")
      txdb <- suppressMessages(makeTxDbFromUCSC(genome=genome))
      print("Making Travis Coordinates ...")
      TravisCoords <- suppressMessages(makeTravisCoordsFromTxDb(txdb))
    } else {
      print("Making Travis Coordinates from provided TranscriptDb Object ...")
      TravisCoords <- makeTravisCoordsFromTxDb(txdb)
    }
  } else {
    print("Using provided Travis Coordinates")
    TravisCoords <- TravisCoordsFromTxDb
  }
  
  # import bed12 file
  noGroup <- length(peak)
  group_names <- names(peak)
  m <- peak
  if (is.null(group_names)) {
    group_names <- paste("item",1:noGroup)
  }
  for (i in 1:noGroup) {
    temp = .countTravisDensity(peak[[i]],TravisCoords)
    temp = cbind(temp,Feature=group_names[i])
    m[[i]] =temp
  }
  ct=.combineListOfDataFrame(m)
  ct[[4]] <- as.character(ct[[4]])
  
  # plot figure
  ct1 <- ct[(ct$category=="mRNA")&(ct$count>0),]
  ct2 <- ct[(ct$category=="lncRNA")&(ct$count>0),]
  
  # add weight
  ct1 <- data.frame(ct1,weight=ct1$count)
  ct2 <- data.frame(ct2,weight=ct2$count)
  featureSet <- as.character(unique(ct1$Feature))
  for (i in 1:length(featureSet)) {
    id <- (ct1$category=="mRNA") & (ct1$Feature==featureSet[i])
    ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
    
    id <- (ct2$category=="mRNA") & (ct2$Feature==featureSet[i])
    ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
    
    id <- (ct1$category=="lncRNA") & (ct1$Feature==featureSet[i])
    ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
    
    id <- (ct2$category=="lncRNA") & (ct2$Feature==featureSet[i])
    ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
  }

  
  
  # save(ct1,ct2,file="TravisPlot.RData")
  p1 <- ggplot(ct1, aes(x=pos, group=Feature, weight=3*weight)) + 
    ggtitle("Distribution on mRNA") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + xlab("") + ylab("Frequency") +
    annotate("pointrange", x = 1, y = -0.3, ymin = -0.5, ymax = -0.1, colour = "black",size=1) + 
    annotate("pointrange", x = 2, y = -0.3, ymin = -0.5, ymax = -0.1, colour = "black",size=1) + 
    annotate("pointrange", x = 0, y = -0.3, ymin = -0.5, ymax = -0.1, colour = "black",size=1) + 
    annotate("pointrange", x = 3, y = -0.3, ymin = -0.5, ymax = -0.1, colour = "black",size=1) + 
    annotate("rect", xmin = 0, xmax = 1, ymin = -0.4, ymax = -0.2, alpha = .2, fill = "blue") +
    annotate("rect", xmin = 1, xmax = 2, ymin = -0.4, ymax = -0.2, alpha = .2, fill = "green") +
    annotate("rect", xmin = 2, xmax = 3, ymin = -0.4, ymax = -0.2, alpha = .2, fill = "red") +
    geom_density(aes(fill=factor(Feature)),alpha=0.5) +
    annotate("text", x = 0.5, y = -0.3, label = "5'UTR") +
    annotate("text", x = 1.5, y = -0.3, label = "CDS") +
    annotate("text", x = 2.5, y = -0.3, label = "3'UTR")
  
  p2 <- ggplot(ct2, aes(x=pos*3, group=Feature, weight=weight*3)) + 
    ggtitle("Distribution on lncRNA") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + xlab("") + ylab("Frequency") +
    annotate("pointrange", x = 0, y = -0.3, ymin = -0.5, ymax = -0.1, colour = "black",size=1) + 
    annotate("pointrange", x = 3, y = -0.3, ymin = -0.5, ymax = -0.1, colour = "black",size=1) + 
    annotate("rect", xmin = 0, xmax = 3, ymin = -0.4, ymax = -0.2, alpha = .2, fill = "green") +
    geom_density(aes(fill=factor(Feature)),alpha=0.5) +
    annotate("text", x = 0.3, y = -0.3, label = "5'End") +
    annotate("text", x = 1.5, y = -0.3, label = "lncRNA") +
    annotate("text", x = 2.7, y = -0.3, label = "3'End") 
  
  if (is.na(saveToPDFprefix)) {
    # return the result
    q <- list(mRNA=p1,lncRNA=p2)
    print("Figures returned as objects ...")
    print("Please check ...")
    return(q)
  }  else {
    f1 <- paste(saveToPDFprefix,"_mRNA.pdf",sep="")
    f2 <- paste(saveToPDFprefix,"_lncRNA.pdf",sep="")
    ggsave(filename=f1,plot=p1,width=6, height=4)
    ggsave(filename=f2,plot=p2,width=6, height=4)
    print(paste("Figures saved into",f1,"and",f2,"...", sep=" "))
  }
  
}


.countTravisDensity <- function(peak,TravisCoords) {
  
  # count overlaps
  n <- countOverlaps(TravisCoords,peak)
  q <- data.frame(mcols(TravisCoords),count=n)
  return(q)
  
}
.combineListOfDataFrame <- function(t){
  txt <- "rbind(t[[1]]"
  for (i in 2:length(t)) {
    txt <- paste(txt,",t[[",i,"]]",sep="")
  }
  txt <- paste(txt,")",sep="")
  c <- parse(text=txt)
  
  newframe <- eval(c)
  
  return(newframe)
}