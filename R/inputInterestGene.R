inputInterestGene <-
function(){
     interestOfGene<-read.table(file.choose(),sep=",",header=FALSE)
	 return(interestOfGene)
}
