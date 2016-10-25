\name{frameCounting}
\alias{frameCounting}

\title{
Counts aligned reads within coding sequence regions by frame and
footprint size, splitting by frame and footprint size.
}
\description{
Ribosome footprint sequencing reads aligning within coding sequence
regions may align in the same frame (relative to start codon) as the
coding sequence, or frame shifted by 1 or 2 frames. This function calls
the number of aligning reads within the coding sequence, split by frame
and footprint size.
}
\usage{
frameCounting(riboDat, cds, lengths = 25:30, offset5p = 0, offset3p = 0)
}
\arguments{
  \item{riboDat}{
    A \code{riboData} object containing the ribosome footprints to be counted. 
}
  \item{cds}{
    A \code{GenomicRanges} object containing the coordinates of the
    coding sequences.
  }
  \item{lengths}{
    Lengths of ribosome footprints to be included in the \code{riboData} object.
  }
  \item{offset5p}{
    An offset to the 5' end of the coding sequence coordinates to
    include ribosome footprints some distance upstream of coding start.}
  \item{offset3p}{
    An offset to the 3' end of the coding sequence coordinates to exclude
    reads too close to coding stop.}
}
\value{
  A \code{riboCoding} object.
}
\author{
Thomas J. Hardcastle
}

\seealso{
  \code{riboData}
}
\examples{

#ribosomal footprint data
datadir <- system.file("extdata", package = "riboSeqR")
ribofiles <- paste(datadir, 
                   "/chlamy236_plus_deNovo_plusOnly_Index", c(17,3,5,7), sep = "")
rnafiles <- paste(datadir, 
                  "/chlamy236_plus_deNovo_plusOnly_Index", c(10,12,14,16), sep = "")

riboDat <- readRibodata(ribofiles, rnafiles, replicates = c("WT", "WT",
"M", "M"))

# CDS coordinates
chlamyFasta <- paste(datadir, "/rsem_chlamy236_deNovo.transcripts.fa", sep = "")
fastaCDS <- findCDS(fastaFile = chlamyFasta, 
                    startCodon = c("ATG"), 
                    stopCodon = c("TAG", "TAA", "TGA"))

# frame calling
fCs <- frameCounting(riboDat, fastaCDS)


}
\keyword{manip}
