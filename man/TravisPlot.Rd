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
TravisPlot(peak, TravisCoordsFromTxDb = NA, txdb = NA, genome = NA, saveToPDFprefix = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak}{
%%     ~~Describe \code{peak} here~~
}
  \item{TravisCoordsFromTxDb}{
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
Jia Meng < jia.meng@xjtlu.edu.cn >
}
\note{
Please let me know if you have any questions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Find the toy data files came with Travis package
bed3=system.file("extdata", "H3K4me3_mm10_1000peaks.bed", package="Travis")
bed12=system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Travis")
bam1=system.file("extdata", "SRR568477_mm10_XIST.bam", package="Travis")
bam2=system.file("extdata", "SRR1768320_mm10_XIST.bam", package="Travis")

# Read 4 different kinds of Genomic Features into R
H3K4me3 <- import.bed(bed3) # bed3 imported as GRanges
m6A <- BED12toGRangesList(bed12) # bed12 imported as GRangesList
SRR568477 <- readGAlignments(bam1) # bam imported as GAlignments
SRR1768320 <- readGAlignments(bam2) # bam imported as GAlignments

# put everything in a list with names
peak <- list(H3K4me3,m6A,SRR568477,SRR1768320) 
names(peak) <- c("H3K4me3", "m6A", "BAM1", "BAM2")

# make TravisPlot
# show the transcriptomic view of 4 sets of genomic features

# OPTION 1
TravisPlot(peak, 
           genome="mm10", 
           saveToPDFprefix = "Toy")
# gene annotation will be downloaded automatically, 
# and two PDF figures with prefix "Toy" will be saved
# in the current working directory

# OPTION 2
# It takes a few minutes to built TravisCoodinatesFromTxDb,
# so we may want to reuse it to save some time.
txdb <- makeTxDbFromUCSC(genome="mm10")
TravisCoordsFromTxDb <- makeTravisCoordsFromTxDb(txdb)
TravisPlot(peak[1:2], 
           TravisCoordsFromTxDb = TravisCoordsFromTxDb, 
           saveToPDFprefix = "Toy2")
TravisPlot(peak[3:4], 
           TravisCoordsFromTxDb = TravisCoordsFromTxDb, 
           saveToPDFprefix = "Toy3")

# OPTION 3
# The function may return the figures which you may resize
tp <- TravisPlot(peak, TravisCoordsFromTxDb = TravisCoordsFromTxDb)
tp[[1]] # Please ignore the warning message
tp[[2]] # Please ignore the warning message
tp[[1]] + ylim(-1,5) # reset the range of y axis
tp[[2]] + ylim(-1,5) # reset the range of y axis
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
