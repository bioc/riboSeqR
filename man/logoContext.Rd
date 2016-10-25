\name{logoContext}
\alias{logoContext}

\title{
Produces a motif of coding start context
}
\description{
The findCDS function reports the immediate context of coding start. This
function uses the 'seqLogo' package to produce a motif of the contexts
reported.
}
\usage{
logoContext(cds, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cds}{
    A GRanges object containing a '$context' value which defines the 7
    bases defining coding start context.
}
  \item{\dots}{
    Other values to be passed to the seqLogo function.
}
}
\value{
  A matrix containing counts of the different bases observed at each
  position of the 7-base context.
}
\author{
Thomas J. Hardcastle
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
