\name{Travis-package}
\alias{Travis-package}
\alias{Travis}
\docType{package}
\title{
Travis
}
\description{
Transcriptomic View of Genomic Features
}
\details{
Transcriptomic View of Genomic Features
}
\author{
Jia Meng
Maintainer: Jia Meng <jia.meng@xjtlu.edu.cn> 
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Install Travis package from Github
library(devtools)
install_github('lzcyzm/Travis')

# load Travis package
library(Travis)

# build travis coordinates from GRangesList or TXDB
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
cds <- cdsBy(txdb)
TravisCoords <- makeTravisCoordsFromGRangesList(cds, noBins=100)
TravisCoordsFromTxDb <- makeTravisCoordsFromTxDb(txdb)

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
# andtwo PDF figures with prefix "Toy" will be saved in the current working directory

# OPTION 2
# It takes a few minutes to built TravisCoodinatesFromTxDb, it can be built separately and reuse
txdb <- makeTxDbFromUCSC(genome="mm10")
TravisCoordsFromTxDb <- makeTravisCoordsFromTxDb(txdb)
TravisPlot(peak, 
           TravisCoordsFromTxDb = TravisCoordsFromTxDb, 
           saveToPDFprefix = "Toy2")

# OPTION 3
# the function can return the two figures which you may resize
tp <- TravisPlot(peak, genome="mm10")
tp[[1]] # Please ignore the warning message initiated by ggplot2
tp[[2]] # Please ignore the warning message initiated by ggplot2
}
