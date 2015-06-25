\name{makeTravisCoordsFromTxDb}
\alias{makeTravisCoordsFromTxDb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Travis Coordinates from TxDb Object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeTravisCoordsFromTxDb(txdb, maximalAmbiguity = 3, minimalComponentLength = 100, minimalNcRNALength = 300, noBins = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txdb}{
%%     ~~Describe \code{txdb} here~~
}
  \item{maximalAmbiguity}{
%%     ~~Describe \code{maximalAmbiguity} here~~
}
  \item{minimalComponentLength}{
%%     ~~Describe \code{minimalComponentLength} here~~
}
  \item{minimalNcRNALength}{
%%     ~~Describe \code{minimalNcRNALength} here~~
}
  \item{noBins}{
%%     ~~Describe \code{noBins} here~~
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
function (txdb, maximalAmbiguity = 3, minimalComponentLength = 100, 
    minimalNcRNALength = 300, noBins = 100) 
{
    require("GenomicFeatures")
    parameter = list()
    parameter$txdb <- txdb
    parameter$maximalAmbiguity <- maximalAmbiguity
    parameter$minimalComponentLength <- minimalComponentLength
    parameter$minimalNcRNALength <- minimalNcRNALength
    parameter$noBins <- noBins
    component <- .extractComponent(parameter)
    print("Building Travis Coordinates. It may take a few minutes ...")
    utr3_bin <- makeTravisCoordsFromGRangesList(component[["utr3"]], 
        parameter$noBins)
    utr5_bin <- makeTravisCoordsFromGRangesList(component[["utr5"]], 
        parameter$noBins)
    cds_bin <- makeTravisCoordsFromGRangesList(component[["cds"]], 
        parameter$noBins)
    ncRNA_bin <- makeTravisCoordsFromGRangesList(component[["ncRNA"]], 
        parameter$noBins)
    print("Travis Coordinates Built ...")
    mcols(utr3_bin) <- data.frame(mcols(utr3_bin), comp = "UTR3", 
        category = "mRNA")
    mcols(utr3_bin)$pos <- mcols(utr3_bin)$pos/3 + 2/3
    mcols(utr5_bin) <- data.frame(mcols(utr5_bin), comp = "UTR5", 
        category = "mRNA")
    mcols(utr5_bin)$pos <- mcols(utr5_bin)$pos/3
    mcols(cds_bin) <- data.frame(mcols(cds_bin), comp = "CDS", 
        category = "mRNA")
    mcols(cds_bin)$pos <- mcols(cds_bin)$pos/3 + 1/3
    mcols(ncRNA_bin) <- data.frame(mcols(ncRNA_bin), comp = "lncRNA", 
        category = "lncRNA")
    TravisCoords <- suppressWarnings(c(utr5_bin, cds_bin, utr3_bin, 
        ncRNA_bin))
    return(TravisCoords)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
