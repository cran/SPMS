\name{txtResultShowing}
\alias{txtResultShowing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To show the text result in R-gui.
}
\description{
The function to show the text results in R-gui.A file of .txt format of all the sub-pathway identification text results is saved in the txtResult folder under the "paths" which is set by the function of
setPath() 
}
\usage{
txtResultShowing(paths, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{name}{
The pathway name that user want to see.
}
}
\value{
It resurns the results.
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

#--txtResultShowing(paths,"hsa")

## The function is currently defined as
function (paths, name) 
{
    if (nchar(name) < 5) {
        n <- sub(" ", "", paste(paths, "/txtResult"))
        setwd(n)
        gR <- readLines(sub(" ", "", paste(paths, "/txtResult/txtResult.txt")))
        gResult <- gsub("\"", "", (gsub("\t", "  ", gR)))
    }
    else {
        n <- sub(" ", "", paste(paths, "/txtResult"))
        setwd(n)
        gR <- readLines(sub(" ", "", paste(paths, "/txtResult/txtResult.txt")))
        gResult1 <- gsub("\"", "", (gsub("\t", "  ", gR)))
        gResult <- data.frame("NA")
        k <- 1
        for (i in 1:length(gResult1)[1]) {
            if (substr(gResult1[i], 1, 8) == name) {
                gResult[k] = gResult1[i]
                k <- k + 1
            }
        }
    }
    return(gResult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ result }
\keyword{ show }% __ONLY ONE__ keyword per line
