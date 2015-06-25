\name{BED12toGRangesList}
\alias{BED12toGRangesList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read BED12 into GRangesList Object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BED12toGRangesList(filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
%%     ~~Describe \code{filepath} here~~
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
function (filepath) 
{
    print("Converting BED12 to GRangesList")
    print("It may take a few minutes")
    a = read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    a = a[, 1:12]
    no_tx = length(a[, 1])
    tx_id = 1:no_tx
    tx_name = paste("line1", 1:no_tx, sep = "")
    tx_chrom = a[, 1]
    tx_strand = a[, 6]
    tx_start = a[, 2] + 1
    tx_end = a[, 3]
    transcripts = data.frame(tx_id, tx_name, tx_chrom, tx_strand, 
        tx_start, tx_end)
    head(transcripts)
    tx_name = tx_name
    gene_id = as.character(a[, 4])
    gene_id[is.na(gene_id)] = "NA"
    gene = data.frame(tx_name, gene_id)
    splicing <- lapply(1:no_tx, .spliceSingleTrans, a = a, tx_start = tx_start)
    splicing <- .combineListOfSplicing(splicing)
    peaks = suppressWarnings(makeTxDb(transcripts = transcripts, 
        splicings = splicing, genes = gene))
    tx <- exonsBy(peaks, "tx", use.names = TRUE)
    mcols(tx) <- a
    return(tx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
