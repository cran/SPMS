\name{deleteAll}
\alias{deleteAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
deleteAll
}
\description{
To delete all folders maked by spftware.
}
\usage{
deleteAll(paths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
}
\references{
None
}
\author{
Xiaomeng Ni
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#--deleteAll(paths)

## The function is currently defined as
function (paths) 
{
    setwd(paths)
    unlink("broad", recursive = TRUE)
    unlink("broad1", recursive = TRUE)
    unlink("coords", recursive = TRUE)
    unlink("depth", recursive = TRUE)
    unlink("depth1", recursive = TRUE)
    unlink("geneAndSymbol", recursive = TRUE)
    unlink("networkData", recursive = TRUE)
    unlink("networkData1", recursive = TRUE)
    unlink("networkData2", recursive = TRUE)
    unlink("nodes", recursive = TRUE)
    unlink("fPathway", recursive = TRUE)
    unlink("picture", recursive = TRUE)
    unlink("temp", recursive = TRUE)
    unlink("txtResult", recursive = TRUE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ delete }
\keyword{ folder }% __ONLY ONE__ keyword per line
