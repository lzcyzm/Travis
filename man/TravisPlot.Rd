\name{TravisPlot}
\alias{TravisPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Transcriptomic View of Genomic Features with Travis Coordinates
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TravisPlot(peak, TravisCoords = NA, txdb = NA, genome = NA, saveToPDFprefix = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak}{
%%     ~~Describe \code{peak} here~~
}
  \item{TravisCoords}{
%%     ~~Describe \code{TravisCoords} here~~
}
  \item{txdb}{
%%     ~~Describe \code{txdb} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{saveToPDFprefix}{
%%     ~~Describe \code{saveToPDFprefix} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (peak, TravisCoords = NA, txdb = NA, genome = NA, saveToPDFprefix = NA) 
{
    require(ggplot2)
    if (is.na(TravisCoords) & is.na(txdb) & is.na(genome)) {
        stop("Most provide one of the three: TravisCoords, txdb or genome")
    }
    if (suppressWarnings(is.na(TravisCoords))) {
        if (suppressWarnings(is.na(txdb))) {
            print("Downloading Transcriptome Information from UCSC ...")
            txdb <- suppressMessages(makeTxDbFromUCSC(genome = genome))
            print("Making Travis Coordinates ...")
            TravisCoords <- suppressMessages(makeTravisCoordsFromTxDb(txdb))
        }
        else {
            print("Making Travis Coordinates from provided TranscriptDb Object ...")
            TravisCoords <- makeTravisCoordsFromTxDb(txdb)
        }
    }
    else {
        print("Using provided Travis Coordinates")
    }
    group_names <- names(peak)
    if (is.null(group_names)) {
        group_names <- paste("item", 1:noGroup)
    }
    noGroup <- length(peak)
    for (i in 1:noGroup) {
        temp = .countTravisDensity(peak[[i]], TravisCoords)
        temp = cbind(temp, Feature = group_names[i])
        peak[[i]] = temp
    }
    ct = .combineListOfDataFrame(peak)
    ct[[5]] <- as.character(ct[[5]])
    ct1 <- ct[ct$category == "mRNA", ]
    ct2 <- ct[ct$category == "lncRNA", ]
    p1 <- ggplot(ct1, aes(x = pos, y = Freq, group = Feature)) + 
        ggtitle("Distribution on mRNA") + theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank()) + xlab("") + ylab("Frequency") + 
        annotate("pointrange", x = 1, y = 0.3, ymin = 0.1, ymax = 0.5, 
            colour = "black", size = 1) + annotate("pointrange", 
        x = 2, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black", 
        size = 1) + annotate("pointrange", x = 0, y = 0.3, ymin = 0.1, 
        ymax = 0.5, colour = "black", size = 1) + annotate("pointrange", 
        x = 3, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black", 
        size = 1) + annotate("rect", xmin = 0, xmax = 1, ymin = 0.15, 
        ymax = 0.4, alpha = 0.2, fill = "blue") + annotate("rect", 
        xmin = 1, xmax = 2, ymin = 0.15, ymax = 0.4, alpha = 0.2, 
        fill = "green") + annotate("rect", xmin = 2, xmax = 3, 
        ymin = 0.15, ymax = 0.4, alpha = 0.2, fill = "red") + 
        geom_line(aes(colour = Feature), size = 2) + annotate("text", 
        x = 0.5, y = 0.3, label = "5'UTR") + annotate("text", 
        x = 1.5, y = 0.3, label = "CDS") + annotate("text", x = 2.5, 
        y = 0.3, label = "3'UTR")
    p2 <- ggplot(ct1, aes(x = pos, y = Freq, group = Feature)) + 
        ggtitle("Distribution on lncRNA") + theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank()) + xlab("") + ylab("Frequency") + 
        annotate("pointrange", x = 0, y = 0.3, ymin = 0.1, ymax = 0.5, 
            colour = "black", size = 1) + annotate("pointrange", 
        x = 3, y = 0.3, ymin = 0.1, ymax = 0.5, colour = "black", 
        size = 1) + annotate("rect", xmin = 0, xmax = 3, ymin = 0.15, 
        ymax = 0.4, alpha = 0.2, fill = "green") + geom_line(aes(colour = Feature), 
        size = 2) + annotate("text", x = 0.3, y = 0.3, label = "5'End") + 
        annotate("text", x = 1.5, y = 0.3, label = "lncRNA") + 
        annotate("text", x = 2.7, y = 0.3, label = "3'End")
    if (is.na(saveToPDFprefix)) {
        q <- list(mRNA = p1, lncRNA = p2)
        print("Figures returned as objects ...")
        print("Please check ...")
        return(q)
    }
    else {
        f1 <- paste(saveToPDFprefix, "_mRNA.pdf", sep = "")
        f2 <- paste(saveToPDFprefix, "_lncRNA.pdf", sep = "")
        ggsave(filename = f1, plot = p1, width = 6, height = 4)
        ggsave(filename = f2, plot = p2, width = 6, height = 4)
        print(paste("Figures saved into", f1, "and", f2, "...", 
            sep = " "))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
