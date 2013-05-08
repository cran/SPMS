\name{getAllSubpathway}
\alias{getAllSubpathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To get all sub-pathways' subgraphs.
}
\description{
Run this function to get all sub-pathways' subgraphs.

Function to save all of graph results of sub-pathway in the pdf file. Files of .pdf format of all sub-pathway graph results are saved in the "picture" folder under the "paths" which is set by the function of
setPath() 

}
\usage{
getAllSubpathway(paths, dORb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{dORb}{
The parameter is convenient to users to choose to run depth or broad funtions.
}
}
\references{
None.
}
\author{
Xiaomeng Ni
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

##--getAllSubpathway(paths,"depth")

## The function is currently defined as
function (paths, dORb) 
{
    print("Please waiting a little time!")
    n1 <- sub(" ", "", paste(paths, "/txtResult"))
    setwd(n1)
    n1 <- sub(" ", "", paste(paths, "/txtResult"))
    setwd(n1)
    gR <- read.table("txtResult.txt", sep = "\t", header = FALSE, 
        stringsAsFactors = FALSE)
    data <- gR[, 9]
    if (dORb == "depth") {
        for (i in 1:length(data)) {
            s1 <- c()
            s2 <- c()
            s3 <- c()
            s1 <- paths
            m <- substr(data[i], (25 + nchar(paths) + 1), (25 + 
                nchar(paths) + 2))
            if (substr(m, 2, 2) != ",") {
                s2 <- as.numeric(m)
            }
            else {
                s2 <- substr(data[i], (25 + nchar(paths) + 1), 
                  (25 + nchar(paths) + 1))
            }
            s3 <- substr(data[i], (28 + nchar(paths) + nchar(s2)), 
                (nchar(data[i]) - 2))
            depth_subPathwayGraph(s1, s2, s3)
        }
    }
    else {
        for (i in 1:length(data)) {
            s1 <- c()
            s2 <- c()
            s3 <- c()
            s1 <- paths
            m <- substr(data[i], (25 + nchar(paths) + nchar(i) + 
                1), (25 + nchar(paths) + nchar(i) + 2))
            if (substr(m, 2, 2) != ",") {
                s2 <- as.numeric(m)
            }
            else {
                s2 <- substr(data[i], (25 + nchar(paths) + nchar(i) + 
                  1), (25 + nchar(paths) + nchar(i) + 1))
            }
            s3 <- substr(data[i], (28 + nchar(paths) + nchar(s2)), 
                (nchar(data[i]) - 2))
            broad_subPathwayGraph(s1, s2, s3)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ subgraph }
\keyword{ all }% __ONLY ONE__ keyword per line