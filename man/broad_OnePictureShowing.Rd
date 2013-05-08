\name{broad_OnePictureShowing}
\alias{broad_OnePictureShowing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To show one picture in R-gui.
}
\description{
Function according to the broad text result to show one result in graphicl way.
}
\usage{
broad_OnePictureShowing(paths, q, pathway)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{q}{
The result's number,according to this parameter to show the result.
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

#--broad_OnePictureShowing(paths,2,'hsa00010')

## The function is currently defined as
function (paths, q, pathway) 
{
    name <- pathway
    pathway <- sub(" ", ".", paste(pathway, "txt"))
    m <- sub(" ", "", paste(paths, "/geneAndSymbol"))
    setwd(m)
    geneAndSymbol <- read.table(pathway, sep = "", header = FALSE)
    n <- sub(" ", "", paste(paths, "/temp/resulttemp"))
    setwd(n)
    gSymbol <- read.delim(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    n1 <- sub(" ", "", paste(paths, "/txtResult"))
    setwd(n1)
    gR <- read.table("txtResult.txt", sep = "\t", header = FALSE)
    gResult <- c()
    for (i in 1:dim(gR)[1]) {
        if (substr(gR[i, 1], 1, nchar(name)) == name) {
            gResult <- rbind(gResult, gR[i, ])
        }
    }
    n2 <- sub(" ", "", paste(paths, "/broad1"))
    setwd(n2)
    broad1 <- read.table(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    n3 <- sub(" ", "", paste(paths, "/broad"))
    setwd(n3)
    broad <- read.table(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    n4 <- sub(" ", "", paste(paths, "/coords"))
    setwd(n4)
    coords <- read.table(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    data(ImageSize, envir = environment())
    ImageSize <- ImageSize
    library(png)
    library(EBImage)
    op <- par(bg = "thistle")
    name1 <- substr(name, (nchar(name) - 4), nchar(name))
    name1 <- paste0("a", name1)
    a <- which(name1 == ImageSize[, 1])
    plot(c(0, as.numeric(as.character(ImageSize[a, 2]))), c(0, 
        as.numeric(as.character(ImageSize[a, 3]))), type = "n", 
        xlab = "", ylab = "", axes = FALSE, ann = FALSE)
    p <- nchar(name)
    species <- c()
    if (p == 8) {
        species <- substr(name, 1, 3)
    }
    else if (p == 9) {
        species <- substr(name, 1, 4)
    }
    p <- sub(" ", "", paste(paths, "/png"))
    setwd(p)
    name1 <- sub(" ", ".", paste(name, "png"))
    image <- readPNG(name1, native = TRUE)
    rasterImage(image, 0, 0, as.numeric(as.character(ImageSize[a, 
        2])), as.numeric(as.character(ImageSize[a, 3])), interpolate = FALSE)
    par(new = TRUE)
    for (i in gSymbol[q, 2]:gSymbol[q, 3]) {
        k <- 0
        k = which(coords[, 1] == broad1[i, gSymbol[q, 1]])
        rect((as.numeric(coords[k, 2]) - 23), abs((as.numeric(coords[k, 
            3]) + 9) - as.numeric(as.character(ImageSize[a, 3]))), 
            (as.numeric(coords[k, 2]) + 23), abs((as.numeric(coords[k, 
                3]) - 9) - as.numeric(as.character(ImageSize[a, 
                3]))), border = "#FF00FF")
    }
    for (i in 10:(9 + as.numeric(as.character(gResult[q, 4])))) {
        for (j in gSymbol[q, 2]:gSymbol[q, 3]) {
            m <- strsplit(broad[j, gSymbol[q, 1]], ", ")[[1]]
            if (is.element(gResult[q, i], m) == TRUE) {
                k <- 0
                k = which(coords[, 1] == broad1[j, gSymbol[q, 
                  1]])
                rect((as.numeric(coords[k, 2]) - 23), abs((as.numeric(coords[k, 
                  3]) + 9) - as.numeric(as.character(ImageSize[a, 
                  3]))), (as.numeric(coords[k, 2]) + 23), abs((as.numeric(coords[k, 
                  3]) - 9) - as.numeric(as.character(ImageSize[a, 
                  3]))), border = "#FF0000")
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ broad }
\keyword{ graphical }% __ONLY ONE__ keyword per line
