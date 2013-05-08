\name{mineSubpathway}
\alias{mineSubpathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use this funtion to mine subpathway.
}
\description{
Function to mine subpathway.
}
\usage{
mineSubpathway(paths, interestOfGene, dORb, name, inG, p_value, broad_path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{interestOfGene}{
The genes of interest that users input.
}
  \item{dORb}{
The parameter is convenient to users to choose to run depth or broad funtions.
}
  \item{name}{
The species name or pathway name.
}
  \item{inG}{
According to the interest genes number(inG) in subpathway to limit the result.
}
  \item{p_value}{
The significance of the p_value.
}
  \item{broad_path}{
The broad_path that user want to limit the result.
}
}
\value{
It returns the final text results.
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

#--mineSubpathway(paths,interestOfGene, "depth","hsa",3,0.01,1)

## The function is currently defined as
function (paths, interestOfGene, dORb, name, inG, p_value, broad_path) 
{
    fPathway <- fPathway(paths, name)
    if (name == "hsa" || substr(name, 1, 3) == "hsa") {
        num <- 39114
    }
    else {
        num <- getNodesNum(paths)
    }
    if (nchar(name) < 5) {
        if (dORb == "depth") {
            data <- depth_subPathwayMinner(paths, num, name, 
                interestOfGene, p_value)
            if (data != "Cannot find subpathway!" && dim(data)[1] > 
                1) {
                data1 <- depth_zhenghe(paths, dORb, inG)
            }
        }
        else {
            data <- broad_subPathwayMinner(paths, num, name, 
                interestOfGene, p_value, broad_path)
            if (data != "Cannot find subpathway!" && dim(data)[1] > 
                1) {
                data1 <- broad_zhenghe(paths, dORb, inG)
            }
        }
    }
    else {
        if (dORb == "depth") {
            data <- depth_subPathwayMinner(paths, num, name, 
                interestOfGene, p_value)
            if (data != "Cannot find subpathway!" && dim(data)[1] > 
                1) {
                data1 <- depth_zhenghe(paths, dORb, inG)
            }
        }
        else {
            data <- broad_subPathwayMinner(paths, num, name, 
                interestOfGene, p_value, broad_path)
            if (data != "Cannot find subpathway!" && data == 
                "Cannot find subpathway!") {
            }
            else {
                data1 <- broad_zhenghe(paths, dORb, inG)
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pathway }
\keyword{ mine }% __ONLY ONE__ keyword per line
