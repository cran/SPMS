deleteAll <-
function(paths){
     setwd(paths)
	 unlink("broad",recursive=TRUE)
     unlink("broad1",recursive=TRUE)
	 unlink("coords",recursive=TRUE)
	 unlink("depth",recursive=TRUE)
	 unlink("depth1",recursive=TRUE)
	 unlink("geneAndSymbol",recursive=TRUE)
	 unlink("networkData",recursive=TRUE)
	 unlink("networkData1",recursive=TRUE)
	 unlink("networkData2",recursive=TRUE)
	 unlink("nodes",recursive=TRUE)
	 unlink("fPathway",recursive=TRUE)
	 unlink("picture",recursive=TRUE)
	 unlink("temp",recursive=TRUE)
	 unlink("txtResult",recursive=TRUE)
	 
}
