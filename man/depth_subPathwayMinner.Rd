\name{depth_subPathwayMinner}
\alias{depth_subPathwayMinner}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
According to the depth traverse datas to mine subpathway.
}
\description{
Function to mine subpathway.
}
\usage{
depth_subPathwayMinner(paths, num, name, interestOfGene, p_value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
The paths is fixed, which is set by setPath function at first.
}
  \item{num}{
The whole nodes number of the species.
}
  \item{name}{
The species name.
}
  \item{interestOfGene}{
The genes that users input to mine subpathway.
}
  \item{p_value}{
The significance of the p_value.
}
}
\value{
It returns the text result that clude the redundant results.
}
\references{
Positional gene enrichment analysis of gene sets for high-resolution identification of overrepresented chromosomal regions
}
\author{
Xiaomeng Ni
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#--depth_subPathwayMinner(paths,num,name,interestOfGene,0.01)

## The function is currently defined as
function (paths, num, name, interestOfGene, p_value) 
{
    if (nchar(name) < 5) {
        k <- nchar(name)
        m <- sub(" ", "", paste(paths, "/fpathway"))
        setwd(m)
        pathway <- paste0(name, ".txt")
        pathwayN <- read.table(pathway, sep = "\t", header = FALSE)
    }
    else {
        pathwayN <- as.matrix(name)
    }
    setwd(paths)
    if (length(which("temp" == dir())) > 0) {
        n <- sub(" ", "", paste(paths, "/temp"))
        setwd(n)
        if (length(dir()) > 0) {
            n <- sub(" ", "", paste(paths, "/temp/resulttemp"))
            setwd(n)
            unlink(dir())
        }
    }
    else {
        dir.create("temp")
    }
    setwd(paths)
    if (length(which("picture" == dir())) > 0) {
        n <- sub(" ", "", paste(paths, "/picture"))
        setwd(n)
        unlink(dir())
    }
    else {
        dir.create("picture")
    }
    x <- sub(" ", "", paste(paths, "/temp"))
    setwd(x)
    if (length(which("unUse1" == dir())) > 0) {
        n <- sub(" ", "", paste(paths, "/temp/unUse1"))
        setwd(n)
        unlink(dir())
    }
    else {
        dir.create("unUse1")
    }
    setwd(x)
    if (length(which("unUse" == dir())) > 0) {
        n <- sub(" ", "", paste(paths, "/temp/unUse"))
        setwd(n)
        unlink(dir())
    }
    else {
        dir.create("unUse")
    }
    print("Please wait a little time!")
    for (z in 1:dim(pathwayN)[1]) {
        pathway <- pathwayN[z, 1]
        setwd(paths)
        a <- "pathwayName.txt"
        pathwayName <- read.table(a, sep = "\t", quote = "", 
            header = FALSE)
        n <- sub(" ", "", paste(paths, "/depth"))
        setwd(n)
        pathway <- sub(" ", ".", paste(pathway, "txt"))
        name1 <- pathway
        name <- substr(pathway, 1, 8)
        pathway <- read.table(pathway, sep = "\t", header = FALSE, 
            stringsAsFactors = FALSE)
        setwd(paths)
        k <- sub(" ", "", paste(paths, "/nodes"))
        setwd(k)
        nodes <- read.table(name1, sep = "\t", header = FALSE, 
            stringsAsFactors = FALSE)
        mark <- c()
        mark1 <- c()
        interestOfGene <- as.vector(as.matrix(interestOfGene))
        for (i in 1:length(interestOfGene)) {
            k <- 0
            k = grep(interestOfGene[i], nodes[, 1])
            if (length(k) > 0) {
                for (j in 1:length(k)) {
                  a <- strsplit(nodes[k[j], 1], ", ")[[1]]
                  if (is.element(interestOfGene[i], a) == TRUE) {
                    mark <- c(mark, nodes[k[j], 1])
                    mark1 <- c(mark1, interestOfGene[i])
                  }
                }
            }
        }
        if (length(mark) < 2) 
            next
        array <- matrix(0, dim(pathway)[1], (dim(pathway)[2])/2)
        array1 <- matrix("0", dim(pathway)[1], (dim(pathway)[2])/2)
        for (i in 1:((dim(pathway)[2])/2)) {
            for (j in 1:length(mark)) {
                k = which(pathway[, (2 * i - 1)] == mark[j])
                if (length(k) > 0) {
                  for (l in 1:length(k)) {
                    array[k[l], i] = 1
                    array1[k[l], i] = mark1[j]
                  }
                }
            }
        }
        g <- num
        q <- length(interestOfGene)
        list1 <- list()
        cl <- c()
        for (p in 1:dim(array)[2]) {
            if (length(diff(which(array[, p] == 1)) - 1) > 0) {
                cl <- c(cl, p)
                ar <- c()
                for (i in 1:length(array[, p])) {
                  if (array[i, p] == 1) {
                    ar <- c(ar, i)
                  }
                }
                col1 <- c()
                for (j in 1:(length(ar) - 1)) {
                  for (l in (j + 1):length(ar)) {
                    col <- c()
                    sum0 <- 0
                    sum1 <- 0
                    col <- cbind(col, ar[j], ar[l])
                    for (k in ar[j]:ar[l]) {
                      if (array[k, p] == 0) {
                        sum0 = sum0 + 1
                      }
                      else {
                        sum1 = sum1 + 1
                      }
                    }
                    col <- cbind(col, sum1, sum0)
                    col <- cbind(col, col[1, 3]/(col[1, 2] - 
                      col[1, 1] + 1))
                    col <- cbind(col, (1 - phyper(col[1, 3], 
                      (col[1, 2] - col[1, 1] + 1), g, q)))
                    col1 <- rbind(col1, col)
                  }
                }
                list1 <- c(list1, list(col1 = col1))
            }
        }
        rule2_data <- list1
        if (length(rule2_data) == 0) 
            next
        if (rule2_data[[1]][1, 6] == "NaN") {
            print("P_value is too small to calculate.")
        }
        for (p in 1:length(cl)) {
            if (dim(rule2_data[[p]])[1] > 1) {
                for (i in 1:(dim(rule2_data[[p]])[1] - 1)) {
                  if (rule2_data[[p]][i, 1] > 0) {
                    for (j in (i + 1):dim(rule2_data[[p]])[1]) {
                      if ((rule2_data[[p]][i, 1] >= rule2_data[[p]][j, 
                        1]) && (rule2_data[[p]][i, 2] <= rule2_data[[p]][j, 
                        2]) && (rule2_data[[p]][i, 4] == rule2_data[[p]][j, 
                        4])) {
                        rule2_data[[p]][i, ] = 0
                      }
                      else if ((rule2_data[[p]][i, 1] <= rule2_data[[p]][j, 
                        1]) && (rule2_data[[p]][i, 2] >= rule2_data[[p]][j, 
                        2]) && (rule2_data[[p]][i, 4] == rule2_data[[p]][j, 
                        4])) {
                        rule2_data[[p]][j, ] = 0
                      }
                    }
                  }
                }
            }
            rule2_data[[p]] <- unique(rule2_data[[p]])
        }
        rule3_data <- rule2_data
        for (p in 1:length(cl)) {
            if (dim(rule3_data[[p]])[1] > 1) {
                for (i in 1:(dim(rule3_data[[p]])[1] - 1)) {
                  if (rule3_data[[p]][i, 1] > 0) {
                    for (j in (i + 1):dim(rule3_data[[p]])[1]) {
                      if ((rule3_data[[p]][i, 1] >= rule3_data[[p]][j, 
                        1]) && (rule3_data[[p]][i, 2] <= rule3_data[[p]][j, 
                        2]) && (rule3_data[[p]][i, 5] < rule3_data[[p]][j, 
                        5])) {
                        rule3_data[[p]][i, ] = 0
                      }
                      else if ((rule3_data[[p]][i, 1] <= rule3_data[[p]][j, 
                        1]) && (rule3_data[[p]][i, 2] >= rule3_data[[p]][j, 
                        2]) && (rule3_data[[p]][i, 5] > rule3_data[[p]][j, 
                        5])) {
                        rule3_data[[p]][j, ] = 0
                      }
                    }
                  }
                }
            }
            rule3_data[[p]] <- unique(rule3_data[[p]])
        }
        rule4_data <- rule3_data
        for (p in 1:length(cl)) {
            if (dim(rule4_data[[p]])[1] > 1) {
                for (i in 1:(dim(rule4_data[[p]])[1] - 1)) {
                  if (rule4_data[[p]][i, 1] > 0) {
                    for (j in (i + 1):dim(rule4_data[[p]])[1]) {
                      if ((rule4_data[[p]][i, 1] >= rule4_data[[p]][j, 
                        1]) && (rule4_data[[p]][i, 2] <= rule4_data[[p]][j, 
                        2]) && (rule4_data[[p]][i, 6] < rule4_data[[p]][j, 
                        6])) {
                        rule4_data[[p]][j, ] = 0
                      }
                      else if ((rule4_data[[p]][i, 1] <= rule4_data[[p]][j, 
                        1]) && (rule4_data[[p]][i, 2] >= rule4_data[[p]][j, 
                        2]) && (rule4_data[[p]][i, 6] > rule4_data[[p]][j, 
                        6])) {
                        rule4_data[[p]][i, ] = 0
                      }
                    }
                  }
                }
            }
            rule4_data[[p]] <- unique(rule4_data[[p]])
        }
        rule5_data <- rule4_data
        for (p in 1:length(cl)) {
            for (i in 1:dim(rule5_data[[p]])[1]) {
                if (q/g * (rule5_data[[p]][i, 2] - rule5_data[[p]][i, 
                  1] + 1) > rule5_data[[p]][i, 3]) {
                  rule5_data[[p]][i, ] = 0
                }
            }
            rule5_data[[p]] <- unique(rule5_data[[p]])
        }
        rule6_data <- rule5_data
        for (p in 1:length(rule6_data)) {
            for (i in 1:dim(rule6_data[[p]])[1]) {
                if (rule6_data[[p]][i, 6] > p_value) {
                  rule6_data[[p]][i, ] = 0
                }
            }
            rule6_data[[p]] <- unique(rule6_data[[p]])
        }
        rule7_data <- rule6_data
        list2 <- list()
        ruledata <- list()
        mark2 <- c()
        for (p in 1:length(rule7_data)) {
            if (max(rule7_data[[p]]) != 0) {
                col2 <- c()
                for (i in 1:dim(rule7_data[[p]])[1]) {
                  if (rule7_data[[p]][i, 1] != 0) {
                    col3 <- c()
                    col3 <- as.character(rule7_data[[p]][i, ])
                    col2 <- rbind(col2, col3)
                    mark2 <- c(mark2, cl[p])
                  }
                }
                list2 <- c(list2, list(col2 = col2))
            }
        }
        mark2 <- unique(mark2)
        ruledata <- list2
        if (length(mark2) == 0) 
            next
        a1 <- matrix(0, dim(pathway)[1], length(mark2))
        a2 <- matrix("0", dim(pathway)[1], length(mark2))
        p1 <- matrix("0", dim(pathway)[1], length(mark2))
        aL <- 1
        for (i in 1:length(mark2)) {
            a1[, aL] <- array[, mark2[i]]
            a2[, aL] <- array1[, mark2[i]]
            p1[, aL] <- pathway[, (2 * mark2[i] - 1)]
            aL <- aL + 1
        }
        arrInG <- a2
        n1 <- sub(" ", "", paste(paths, "/depth1"))
        setwd(n1)
        depth1 <- read.table(name1, header = FALSE, sep = "\t", 
            na.strings = "NA", stringsAsFactors = FALSE)
        for (p in 1:length(ruledata)) {
            for (i in 1:dim(ruledata[[p]])[1]) {
                for (j in ruledata[[p]][i, 1]:(as.numeric(ruledata[[p]][i, 
                  2]) - 1)) {
                  k <- 0
                  col <- c()
                  k = which(depth1[1, ] == depth1[j, (2 * mark2[p] - 
                    1)])
                  for (l in (j + 1):ruledata[[p]][i, 2]) {
                    r <- 0
                    r = which(depth1[, k] == depth1[l, (2 * mark2[p] - 
                      1)])
                    col <- c(col, depth1[r, (k + 1)])
                  }
                  if (min(col) > 1) {
                    ruledata[[p]][i, ] = 0
                    break
                  }
                }
            }
        }
        list3 <- list()
        for (p in 1:length(ruledata)) {
            for (i in 1:dim(ruledata[[p]])[1]) {
                if (ruledata[[p]][i, 1] != 0) {
                  col4 <- c()
                  col4 <- cbind(col4, (2 * mark2[p] - 1))
                  col4 <- cbind(col4, ruledata[[p]][i, 1], ruledata[[p]][i, 
                    2])
                  for (j in ruledata[[p]][i, 1]:ruledata[[p]][i, 
                    2]) {
                    col4 <- cbind(col4, as.character(p1[j, p]))
                  }
                  list3 <- c(list3, list(col4 = col4))
                }
            }
        }
        if (length(list3) == 0) 
            next
        data1 <- list3
        datan <- data.frame()
        for (p in 1:length(data1)) {
            for (i in 1:length(data1[[p]])) {
                datan[p, i] = data1[[p]][i]
            }
        }
        setwd(paths)
        m <- sub(" ", "", paste(paths, "/temp/unUse1"))
        setwd(m)
        write.table(datan, name1, quote = T, sep = "\t", col.names = FALSE, 
            append = TRUE, row.names = FALSE, na = "NA")
        number <- pathwayName[, 1]
        no <- sub("a", "", number)
        list4 <- list()
        sum <- 0
        for (p in 1:length(ruledata)) {
            for (i in 1:(length(ruledata[[p]])/6)) {
                if (ruledata[[p]][i, 1] != 0) {
                  col5 <- c()
                  k <- 0
                  col5 <- cbind(col5, name, i + sum)
                  k = which(substr(name, 4, 8) == no)
                  col5 <- cbind(col5, as.character(pathwayName[k, 
                    2]))
                  col5 <- cbind(col5, ruledata[[p]][i, 3], (as.numeric(ruledata[[p]][i, 
                    2]) - as.numeric(ruledata[[p]][i, 1]) + 1), 
                    length(interestOfGene), g, ruledata[[p]][i, 
                      6])
                  subFunc <- gsub(" ", "", paste("depth_subPathwayGraph(", 
                    "'", paths, "'", ",", (i + sum), ",", "'", 
                    name, "'", ")"))
                  col5 <- cbind(col5, subFunc)
                  for (j in ruledata[[p]][i, 1]:ruledata[[p]][i, 
                    2]) {
                    if (arrInG[j, p] != 0) {
                      col5 <- cbind(col5, as.character(arrInG[j, 
                        p]))
                    }
                  }
                  list4 <- c(list4, list(col5 = col5))
                }
            }
            sum = sum + (length(ruledata[[p]])/6)
        }
        data5 <- data.frame()
        for (p in 1:length(list4)) {
            for (i in 1:length(list4[[p]])) {
                data5[p, i] = list4[[p]][i]
            }
        }
        setwd(paths)
        k <- sub(" ", "", paste(paths, "/temp/unUse"))
        setwd(k)
        write.table(data5, name1, quote = T, sep = "\t", col.names = FALSE, 
            append = TRUE, row.names = FALSE, na = "NA")
    }
    if (exists("data5") && length(data5) > 0) {
        return(data5)
    }
    else {
        print("Cannot find subpathway!")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mine }
\keyword{ pathway }% __ONLY ONE__ keyword per line
