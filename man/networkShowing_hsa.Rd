\name{networkShowing_hsa}
\alias{networkShowing_hsa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To show the network in R-gui.
}
\description{
Function to show the network in R-gui.
}
\usage{
networkShowing_hsa(paths, pathway)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{pathway}{
The species name of user wants to see.
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

#--networkShowing_hsa(paths,"hsa00052")

## The function is currently defined as
function (paths, pathway) 
{
    library(igraph)
    n <- sub(" ", "", paste(paths, "/networkData2"))
    setwd(n)
    name <- gsub(" ", "", paste(pathway, ".txt"))
    data <- read.table(name, sep = "\t", header = FALSE)
    d = data.frame(p1 = data[, 1], p2 = data[, 2])
    g = graph.data.frame(d, directed = FALSE)
    tkplot(g, vertex.label = V(g)$name, vertex.size = 20)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ network }
\keyword{ show }% __ONLY ONE__ keyword per line
