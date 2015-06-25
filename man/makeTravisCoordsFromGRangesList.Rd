\name{makeTravisCoordsFromGRangesList}
\alias{makeTravisCoordsFromGRangesList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Travis Coordinates from GRangesList Object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeTravisCoordsFromGRangesList(comp, noBins)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comp}{
%%     ~~Describe \code{comp} here~~
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
function (comp, noBins) 
{
    require("GenomicFeatures")
    tx_length <- as.numeric(sum(width(comp)))
    checkpoints_interval <- (tx_length - 1)/(noBins - 1)
    tx_name <- names(comp)
    granges <- unlist(comp)
    tx <- granges[tx_name]
    strand <- as.character(as.data.frame(strand(tx))[[1]])
    chr <- as.character(as.data.frame(seqnames(tx))[[1]])
    t <- lapply(X = 1:noBins, FUN = .makeTravisCoordsForSingleIndex, 
        comp = comp, noBins = noBins, checkpoints_interval = checkpoints_interval, 
        strand = strand, chr = chr, tx_name = tx_name)
    TravisCoords <- .combineListOfGRanges(t)
    return(TravisCoords)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
