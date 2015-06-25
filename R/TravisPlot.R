

TravisPlot <- function(peak,TravisCoords=NA,txdb=NA,genome=NA,saveToPDFprefix=NA){
  
  # make sure the travis coordinates are available
  if (is.na(TravisCoords)&is.na(txdb)&is.na(genome)) {
    stop("Most provide one of the three: TravisCoords, txdb or genome")
  } 
  
  if ( suppressWarnings(is.na(TravisCoords)) ) {
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
  }
  
  # import bed12 file
  noGroup <- length(peak)
  group_names <- names(peak)
  if (is.null(group_names)) {
    group_names <- paste("item",1:noGroup)
  }
  for (i in 1:noGroup) {
    temp = .countTravisDensity(peak[[i]],TravisCoords)
    temp = cbind(temp,Feature=group_names[i])
    peak[[i]] =temp
  }
  ct=.combineListOfDataFrame(peak)
  ct[[5]] <- as.character(ct[[5]])
  
  # plot figure
  ct1 <- ct[ct$category=="mRNA",]
  ct2 <- ct[ct$category=="lncRNA",]
  
  p1 <- ggplot(ct1, aes(x=pos, y=Freq, group=Feature)) + 
    ggtitle("Distribution on mRNA") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + xlab("") + ylab("Frequency") +
    annotate("pointrange", x = 1, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black",size=1) + 
    annotate("pointrange", x = 2, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black",size=1) + 
    annotate("pointrange", x = 0, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black",size=1) + 
    annotate("pointrange", x = 3, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black",size=1) + 
    annotate("rect", xmin = 0, xmax = 1, ymin = 0.15, ymax = 0.4, alpha = .2, fill = "blue") +
    annotate("rect", xmin = 1, xmax = 2, ymin = 0.15, ymax = 0.4, alpha = .2, fill = "green") +
    annotate("rect", xmin = 2, xmax = 3, ymin = 0.15, ymax = 0.4, alpha = .2, fill = "red") +
    geom_line(aes(colour = Feature),size=2) +
    annotate("text", x = 0.5, y = 0.3, label = "5'UTR") +
    annotate("text", x = 1.5, y = 0.3, label = "CDS") +
    annotate("text", x = 2.5, y = 0.3, label = "3'UTR")
  
  p2 <- ggplot(ct1, aes(x=pos, y=Freq, group=Feature)) + 
    ggtitle("Distribution on lncRNA") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + xlab("") + ylab("Frequency") +
    annotate("pointrange", x = 0, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black",size=1) + 
    annotate("pointrange", x = 3, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black",size=1) + 
    annotate("rect", xmin = 0, xmax = 3, ymin = 0.15, ymax = 0.4, alpha = .2, fill = "green") +
    geom_line(aes(colour = Feature),size=2) +
    annotate("text", x = 0.3, y = 0.3, label = "5'End") +
    annotate("text", x = 1.5, y = 0.3, label = "lncRNA") +
    annotate("text", x = 2.7, y = 0.3, label = "3'End") 
  
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
  
  # divide by mRNA and lncRNA
  mrna <- q[(q$category=="mRNA")&(q$count>0),]
  lncrna <- q[(q$category=="lncRNA")&(q$count>0),]
  
  # process mRNA
  mrna$comp <- as.character(mrna$comp)
  result <- table(mrna[,2:3])
  mrna <- as.data.frame(result)
  mrna$Freq <- mrna$Freq/mean(mrna$Freq)/3
  mrna <- data.frame(mrna,category="mRNA")
  mrna <- mrna[(mrna$Freq>0.0000000001),]
  
  
  # process lncRNA
  lncrna$comp <- as.character(lncrna$comp)
  result <- table(lncrna[,c(2,3)])
  lncrna <- as.data.frame(result)
  lncrna$Freq <- lncrna$Freq/mean(lncrna$Freq)
  lncrna <- data.frame(lncrna,category="lncRNA")
  
  
  # output
  out <- rbind(mrna,lncrna)
  out$pos <- as.numeric(as.character(out$pos))
  out$comp <- factor(as.character(out$comp),levels=c("UTR5","CDS","UTR3","lncRNA"))
  out$pos <- out$pos*3
  
  # return
  return(out)
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