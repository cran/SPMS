\name{fPathway}
\alias{fPathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain  the number of pathway in one species.
}
\description{
This function returns a pathway details of one species and will be used later.
}
\usage{
fPathway(paths, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{name}{
The name of one species. 
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

#--fPathway(paths,"hsa")

## The function is currently defined as
function (paths, name) 
{
    setwd(paths)
    if (length(which("fPathway" == dir())) > 0) {
        n <- sub(" ", "", paste(paths, "/fPathway"))
        setwd(n)
        unlink(dir())
    }
    else {
        dir.create("fPathway")
    }
    sum <- 0
    m <- sub(" ", "", paste(paths, "/depth"))
    setwd(m)
    k <- nchar(name)
    pathwayN <- c()
    for (i in 1:length(dir())) {
        m <- sub(" ", "", paste(paths, "/depth"))
        setwd(m)
        if (substr(dir()[i], 1, k) == name) {
            pathwayN <- rbind(pathwayN, substr(dir()[i], 1, (k + 
                5)))
        }
    }
    n <- sub(" ", "", paste(paths, "/fPathway"))
    setwd(n)
    name <- sub(" ", "", paste(name, ".txt"))
    write.table(pathwayN, name, col.names = FALSE, row.names = FALSE)
    sum = sum + 1
    return(sum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pathway }
\keyword{ number }% __ONLY ONE__ keyword per line
