getAllSubpathway <-
function(paths,dORb){
     print("Please waiting a little time!")
	 n1<-sub(" ","",paste(paths,"/txtResult"))
	 setwd(n1)
	 n1<-sub(" ","",paste(paths,"/txtResult"))
	 setwd(n1)
	 gR<-read.table("txtResult.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
     data<-gR[,9]
	 if(dORb=="depth"){
	     for(i in 1:length(data)){
		     s1<-c()
		     s2<-c()
		     s3<-c()
		     s1<-paths
			 m<-substr(data[i],(25+nchar(paths)+1),(25+nchar(paths)+2))
			 if(substr(m,2,2)!=","){
			     s2<-as.numeric(m)
			 }else{
			     s2<-substr(data[i],(25+nchar(paths)+1),(25+nchar(paths)+1))
			 }
			 s3<-substr(data[i],(28+nchar(paths)+nchar(s2)),(nchar(data[i])-2))
			 depth_subPathwayGraph(s1,s2,s3)
		 }
	 }else{
	     for(i in 1:length(data)){
		     s1<-c()
		     s2<-c()
		     s3<-c()
		     s1<-paths
			 m<-substr(data[i],(25+nchar(paths)+1),(25+nchar(paths)+2))
			 if(substr(m,2,2)!=","){
			     s2<-as.numeric(m)
			 }else{
			     s2<-substr(data[i],(25+nchar(paths)+1),(25+nchar(paths)+1))
			 }
			 s3<-substr(data[i],(28+nchar(paths)+nchar(s2)),(nchar(data[i])-2))
			 broad_subPathwayGraph(s1,s2,s3)
		 }
	 }
	  
}
