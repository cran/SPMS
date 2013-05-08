\name{setPath}
\alias{setPath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To set path which used to save all datas and results.
}
\description{
Function to set the path that the whole data and result will store.
}
\usage{
setPath()
}
\value{
It returns the paths.
}
\references{
Nnone
}
\author{
Xiaomeng Ni
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## This function must return a path and named it as paths
#--paths<-setPath()

## The function is currently defined as
setPath<-function () 
{
    paths <- scan(what = "character")
    return(paths)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ setpath }
\keyword{ pathway }% __ONLY ONE__ keyword per line
