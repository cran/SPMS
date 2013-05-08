\name{getNodesNum}
\alias{getNodesNum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The whole nodes number of the species.
}
\description{
To get the whole nodes number of the species.
}
\usage{
getNodesNum(paths)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
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

#--num<-getNodesNum(paths)

## The function is currently defined as
function (paths) 
{
    library(XML)
    n1 <- paste0(paths, "/xml")
    setwd(n1)
    num <- 0
    for (p in 1:length(dir())) {
        n <- paste0(paths, "/xml")
        setwd(n)
        pathway <- dir()[p]
        top <- xmlRoot(xmlTreeParse(pathway))
        KOList1 <- character()
        KOListIndex1 <- 0
        gene_symbol <- character()
        for (i in 1:length(top)) {
            if (xmlName(top[[i]]) == "entry") {
                if (!is.na(xmlAttrs(top[[i]])["type"])) {
                  if (!is.na(xmlAttrs(top[[i]])["name"])) {
                    if (xmlAttrs(top[[i]])["type"] == "gene") {
                      tmp <- unlist(xmlChildren(top[[i]]))
                      if (length(tmp) > 0) {
                        for (k in 1:length(tmp)) {
                          KOListIndex1 <- KOListIndex1 + 1
                          KOList1[KOListIndex1] <- tmp[k]
                        }
                      }
                    }
                  }
                }
            }
        }
        gene_symbol <- KOList1
        gene_AllSymbol <- gene_symbol
        gene_symbol <- c()
        for (i in 1:length(gene_AllSymbol)) {
            if (gene_AllSymbol[i] == "graphics") {
                gene_symbol <- c(gene_symbol, gene_AllSymbol[i + 
                  1])
            }
        }
        num <- num + length(gene_symbol)
    }
    return(num)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nodes }
\keyword{ number }% __ONLY ONE__ keyword per line
