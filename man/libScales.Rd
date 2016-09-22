\name{libScales}
\alias{libScales}
\title{
Generates library scaling factors from ribosome footprint (and
associated RNA) data.
}
\description{
This function applies the library scaling factor estimation methods
implemented in the baySeq package to ribosome footprint data and (if
present) the associated RNA-Seq data. The library scaling factors can
then be used in downstream analyses. 
}
\usage{
libScales(rC, riboDat, lengths, frames, method = "edgeR")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rC}{
    A \code{\link{riboCoding}} object containing the coordinates of the
    coding sequences and the number of ribosomal footprint reads of
    various classes to be found in each.
  }
  \item{riboDat}{
    A \code{\link{riboData}} object containing the ribosome footprint
    and RNA-seq alignments.
  }
  \item{lengths}{
    Lengths of ribosome footprints to inform count data.
  }
  \item{frames}{
    Frames of ribosome footprints (relative to coding start site). If
    omitted, all frames are used. 
}
  \item{method}{
    Method to be used for scaling factor estimation. See \code{\link[baySeq]{getLibsizes}}
}
}
\value{
A list object containing
\item{riboLS}{Ribosome footprint library scaling factors.}
\item{rnaLS}{RNA-Seq library scaling factors.}
}
\author{
Thomas J. Hardcastle
}
\keyword{manip}
