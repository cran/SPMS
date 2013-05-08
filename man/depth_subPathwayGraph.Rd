\name{depth_subPathwayGraph}
\alias{depth_subPathwayGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
To save the graphicl result in pdf format.
}
\description{
Function to show the text result of depth in picture.
}
\usage{
depth_subPathwayGraph(paths, q, pathway)
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

##-depth_subPathwayGraph(paths,1,'hsa00052')

## The function is currently defined as
function (paths, q, pathway) 
{
    q <- as.numeric(q)
    name <- pathway
    pathway <- sub(" ", ".", paste(pathway, "txt"))
    m <- sub(" ", "", paste(paths, "/geneAndSymbol"))
    setwd(m)
    geneAndSymbol <- read.table(pathway, sep = "", header = FALSE, 
        stringsAsFactors = FALSE)
    n <- sub(" ", "", paste(paths, "/temp/resulttemp"))
    setwd(n)
    gSymbol <- read.delim(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    n1 <- sub(" ", "", paste(paths, "/txtResult"))
    setwd(n1)
    gR <- read.table("txtResult.txt", sep = "\t", header = FALSE, 
        stringsAsFactors = FALSE)
    gResult <- c()
    for (i in 1:dim(gR)[1]) {
        if (substr(gR[i, 1], 1, nchar(name)) == name) {
            gResult <- rbind(gResult, gR[i, ])
        }
    }
    n2 <- sub(" ", "", paste(paths, "/depth1"))
    setwd(n2)
    depth1 <- read.table(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    n3 <- sub(" ", "", paste(paths, "/depth"))
    setwd(n3)
    depth <- read.table(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    n4 <- sub(" ", "", paste(paths, "/coords1"))
    setwd(n4)
    coords <- read.table(pathway, header = FALSE, sep = "\t", 
        na.strings = "NA", stringsAsFactors = FALSE)
    setwd(paths)
    size <- read.delim("ImageSize.txt", header = FALSE, na.strings = "NA", 
        sep = "\t")
    library(png)
    o <- sub(" ", "", paste(paths, "/picture"))
    setwd(o)
    pdffile <- paste0(gResult[q, 1], ".pdf")
    pdf(file = pdffile, bg = "transparent")
    library(EBImage)
    op <- par(bg = "thistle")
    a <- which(name == size[, 1])
    plot(c(0, as.numeric(as.character(size[a, 2]))), c(0, as.numeric(as.character(size[a, 
        3]))), type = "n", xlab = "", ylab = "", axes = FALSE, 
        ann = FALSE)
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
    rasterImage(image, 0, 0, as.numeric(as.character(size[a, 
        2])), as.numeric(as.character(size[a, 3])), interpolate = FALSE)
    par(new = TRUE)
    for (i in gSymbol[q, 2]:gSymbol[q, 3]) {
        k <- 0
        k = which(coords[, 1] == depth1[i, gSymbol[q, 1]])
        rect((as.numeric(coords[k, 2]) - 23), abs((as.numeric(coords[k, 
            3]) + 9) - as.numeric(as.character(size[a, 3]))), 
            (as.numeric(coords[k, 2]) + 23), abs((as.numeric(coords[k, 
                3]) - 9) - as.numeric(as.character(size[a, 3]))), 
            border = "#FF00FF")
    }
    for (i in 10:(9 + gResult[q, 4])) {
        for (j in gSymbol[q, 2]:gSymbol[q, 3]) {
            m <- strsplit(depth[j, gSymbol[q, 1]], ", ")[[1]]
            if (is.element(gResult[q, i], m) == TRUE) {
                k <- 0
                k = which(coords[, 1] == depth1[j, gSymbol[q, 
                  1]])
                rect((as.numeric(coords[k, 2]) - 23), abs((as.numeric(coords[k, 
                  3]) + 9) - as.numeric(as.character(size[a, 
                  3]))), (as.numeric(coords[k, 2]) + 23), abs((as.numeric(coords[k, 
                  3]) - 9) - as.numeric(as.character(size[a, 
                  3]))), border = "#FF0000")
            }
        }
    }
    dev.off()
    getwd()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }
\keyword{ pathway }% __ONLY ONE__ keyword per line
