fPathway <-
function(paths,name){
     setwd(paths)
	 if(length(which("fPathway"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/fPathway")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("fPathway")
	 }
     sum<-0	 
	     m<-sub(" ","",paste(paths,"/depth"))
	     setwd(m)
	     k<-nchar(name)
		 pathwayN<-c()
	     for(i in 1:length(dir())){
		     m<-sub(" ","",paste(paths,"/depth"))
	         setwd(m) 
	         if(substr(dir()[i],1,k)==name){
		         pathwayN<-rbind(pathwayN,substr(dir()[i],1,(k+5)))
		     }
	     }
		 n<-sub(" ","",paste(paths,"/fPathway"))
	     setwd(n) 
		 name<-sub(" ","",paste(name,".txt"))
		 write.table(pathwayN,name,col.names=FALSE,row.names=FALSE)
		 sum=sum+1
	 return(sum) 
}
