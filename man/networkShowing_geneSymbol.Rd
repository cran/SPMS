\name{networkShowing_geneSymbol}
\alias{networkShowing_geneSymbol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To show the network in R-gui.
}
\description{
Function to show the network in R-gui.
}
\usage{
networkShowing_geneSymbol(paths, pathway)
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

##--networkShowing_geneSymbol(paths,"hsa00052")

## The function is currently defined as
function (paths, pathway) 
{
    library(igraph)
    n <- sub(" ", "", paste(paths, "/networkData1"))
    setwd(n)
    name <- gsub(" ", "", paste(pathway, ".txt"))
    data2 <- read.table(name, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    n <- sub(" ", "", paste(paths, "/networkData"))
    setwd(n)
    data3 <- read.table(name, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    data <- c(data2[, 1], data2[, 2])
    data1 <- c(data3[, 1], data3[, 2])
    for (i in 1:length(data)) {
        k = which(data[i] == data)
        if (length(k) > 1) {
            for (j in 1:(length(k) - 1)) {
                for (l in (j + 1):length(k)) {
                  if (data1[k[j]] != data1[k[l]]) {
                    data[k[l]] = paste(data[k[l]], "")
                  }
                }
                break
            }
        }
    }
    data <- cbind(data[1:(length(data)/2)], data[(length(data)/2 + 
        1):length(data)])
    d = data.frame(p1 = data[, 1], p2 = data[, 2])
    g = graph.data.frame(d, directed = FALSE)
    tkplot(g, vertex.label = V(g)$name, vertex.size = 20)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ network }
\keyword{ show }% __ONLY ONE__ keyword per line
