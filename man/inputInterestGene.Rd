\name{inputInterestGene}
\alias{inputInterestGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Input interestGene.
}
\description{
Users can input to interest genes according to the path. The genes input by user can be perpared in one file and each line represents one gene. The gene IDs identified  by this function are the same as the "graphics name" in xml files and they may be the gene symbol ID or ORF nameID ,etc. 
}
\usage{
inputInterestGene()
}
%- maybe also 'usage' for other objects documented here.
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

##This function's returnd value must named as interestOfGene.
#--interestOfGene<-inputInterestGene()

## The function is currently defined as
function () 
{
    interestOfGene <- read.table(file.choose(), sep = ",", header = FALSE)
    return(interestOfGene)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interest }
\keyword{ gene }% __ONLY ONE__ keyword per line
