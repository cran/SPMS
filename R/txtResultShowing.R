txtResultShowing <-
function(paths,name){
     if(nchar(name)<5){
	     n<-sub(" ","",paste(paths,"/txtResult"))
	     setwd(n)
	     gR<-readLines(sub(" ","",paste(paths,"/txtResult/txtResult.txt")))
         gResult<-gsub("\"","",(gsub("\t","  ",gR)))
	 }else{
	     n<-sub(" ","",paste(paths,"/txtResult"))
	     setwd(n)
	     gR<-readLines(sub(" ","",paste(paths,"/txtResult/txtResult.txt")))
         gResult1<-gsub("\"","",(gsub("\t","  ",gR)))
		 gResult<-data.frame("NA")
		 k<-1
	     for(i in 1:length(gResult1)[1]){
             if(substr(gResult1[i],1,8)==name){
			     gResult[k]=gResult1[i]
				 k<-k+1
			 }
	     }
	 }
	 
	 return(gResult)
}
