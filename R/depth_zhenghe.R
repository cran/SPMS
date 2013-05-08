depth_zhenghe <-
function(paths,dORb,inG){
   w<-sub(" ","",paste(paths,"/temp"))
   setwd(w)
   if(length(which("resulttemp"==dir()))>0){
		 n<-sub(" ","",paste(w,"/resulttemp")) 
	     setwd(n)
		 unlink(dir())
   }else{
		 dir.create("resulttemp")
   }
   setwd(paths)
   if(length(which("txtResult"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/txtResult")) 
	     setwd(n)
		 unlink(dir())
   }else{
		 dir.create("txtResult")
   }
   num<-1
   a<-sub(" ","",paste(paths,"/temp/unUse"))
   setwd(a)
   for(z in 1:length(dir(a))){
     if(length(dir(a))==0)break
     name<-substr(dir(a)[1],1,8)
	 pathway<-dir(a)[1]
	 a<-sub(" ","",paste(paths,"/temp/unUse1"))
	 setwd(a)	 
	 result11<-read.delim(pathway,sep="\t",header=FALSE,quote="",na.strings="NA",stringsAsFactors=FALSE)
	 unlink(pathway)
	 b<-sub(" ","",paste(paths,"/temp/unUse"))
	 setwd(b)	 
	 result10<-read.table(pathway,sep="\t",header=FALSE,na.strings="NA",stringsAsFactors=FALSE)
	 unlink(pathway)
	
     result1<-result11
	 result<-result10
	 if(dim(result)[1]>1){
	  for(i in 1:(dim(result)[1]-1)){		     
	     for(j in (i+1):dim(result)[1]){
             if(result[i,4]==result[j,4]&&result[i,5]==result[j,5]){
			     n<-c()
		         m<-c()
			     for(l in 4:(4+result[i,5]-1)){				     
				     n<-c(n,as.character(result1[i,l]))					 
			     }
		         for(p in 4:(4+result[j,5]-1)){
				     m<-c(m,as.character(result1[j,p]))
			     }
			     if(length(pmatch(n,m))==result[i,5]&&length(pmatch(n,m))==result[j,5]&&(is.na(sum(pmatch(n,m)))==FALSE)){
				     result[j,]=0
				     result1[j,]=0
			     }
			 }		 			 			 
		 }
	  }
	 }
	 result<-unique(result)
	 result1<-unique(result1)
	 
	 if(dim(result)[1]>1){
	  for(i in 1:(dim(result)[1]-1)){		     
	     for(j in (i+1):dim(result)[1]){
			 if(result[i,4]==result[j,4]&&result[i,5]<=result[j,5]){
			     n<-c()
				 m<-c()
			     for(l in 10:(10+result[i,4]-1)){				     
				     n<-c(n,as.character(result[i,l]))					 
				 }
				 for(p in 10:(10+result[j,4]-1)){
				     m<-c(m,as.character(result[j,p]))
				 }
				 if(length(pmatch(n,m))==result[i,4]&&(is.na(sum(pmatch(n,m)))==FALSE)){
				     result[j,]=0
					 result1[j,]=0
				 }
			 }
             if(result[i,4]==result[j,4]&&result[i,5]>=result[j,5]){
			     n<-c()
				 m<-c()
			     for(l in 10:(10+result[i,4]-1)){				     
				     n<-c(n,as.character(result[i,l]))					 
				 }
				 for(p in 10:(10+result[j,4]-1)){
				     m<-c(m,as.character(result[j,p]))
				 }
				 if(length(pmatch(m,n))==result[j,4]&&(is.na(sum(pmatch(m,n)))==FALSE)){
				     result[i,]=0
					 result1[i,]=0
				 }
			 }			 
		 }			 
	  }
	 } 
	 result<-unique(result)
	 result1<-unique(result1)
	 
	 if(dim(result)[1]>1){
	  for(i in 1:(dim(result)[1]-1)){		     
	     for(j in (i+1):dim(result)[1]){
			 if(result[i,4]<result[j,4]&&((result[i,5]-result[i,4])==(result[j,5]-result[j,4]))){
			     n<-c()
				 m<-c()
			     for(l in 4:(4+result[i,5]-1)){				     
				     n<-c(n,as.character(result1[i,l]))					 
				 }
				 for(p in 4:(4+result[j,5]-1)){
				     m<-c(m,as.character(result1[j,p]))
				 }
				 if(length(pmatch(n,m))==result[i,5]){
				     result[i,]=0
					 result1[i,]=0
				 }
			 }
             if(result[i,4]>result[j,4]&&((result[i,5]-result[i,4])==(result[j,5]-result[j,4]))){
			     n<-c()
				 m<-c()
			     for(l in 4:(4+result[i,5]-1)){				     
				     n<-c(n,as.character(result1[i,l]))					 
				 }
				 for(p in 4:(4+result[j,5]-1)){
				     m<-c(m,as.character(result1[j,p]))
				 }
				 if(length(pmatch(m,n))==result[j,5]){
				     result[j,]=0
					 result1[j,]=0
				 }
			 }			 
		 }			 
	  }
	 }
	 result<-unique(result)
	 result1<-unique(result1)
	 
	 for(i in 1:dim(result)[1]){
	     if(result[i,4]<inG){
		     result[i,]=0
			 result1[i,]=0
		 }
	 }
     result<-unique(result)
	 result1<-unique(result1)
	 if(max(result[,4])==0)next
	
	 for(i in 1:dim(result)[1]){
	     mark1<-c()
		 for(j in 1:result[i,4]){
	         mark1<-c(mark1,result[i,(9+j)])
	     }
		 if(length(unique(mark1))==1){
		     result[i,]=0
			 result1[i,]=0
		 }else{
		     mark2<-c()
			 mark2<-unique(mark1)
			 for(j in 1:length(mark2)){
			     k<-0
				 k=which(mark1==mark2[j])
				 if(length(k)/length(mark1)>=0.5){
				     result[i,]=0
			         result1[i,]=0
					 break
				 }
			 }
		 }
	 }
	 result<-unique(result)
	 result1<-unique(result1)
	 
	 for(i in 1:dim(result)[1]){
	     if(result[i,1]==0){
		     result<-result[-i,]
			 result1<-result1[-i,]
			 break
		 }
	 }
	 if(dim(result)[1]==0)next
	 
	 for(i in 1:dim(result)[1]){
	     result[i,2]=num
         result[i,9]<-as.character(gsub(" ","",paste("depth_subPathwayGraph(","\'",paths,"\'",",",i,",","\'",name,"\'",")")))		 
		 num<-num+1
	 }
     c<-sub(" ","",paste(paths,"/txtResult"))
	 setwd(c)
	 write.table(result,"txtResult.txt",sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,quote=TRUE)
	 d<-sub(" ","",paste(paths,"/temp/resulttemp"))
	 setwd(d)
	 write.table(result1,pathway,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,na="NA",quote=FALSE)
  }
     
     #====================================================================	 
	 n1<-sub(" ","",paste(paths,"/txtResult"))
	 setwd(n1)
   if(length(dir(n1))>0){
	 gR<-scan("txtResult.txt",what=character(),sep="\t",na.strings="NA")
	 unlink("txtResult.txt")
     k<-c()
     gResult1<-data.frame()
	 k<-c(k,grep("hsa0",gR))
     if(length(k)>2){
         for(i in 1:(length(k)/2-1)){
             l<-1
        	 gResult2<-matrix("0",1,50)
             for(j in k[2*i-1]:(k[2*i+1]-1)){
        	     gResult2[1,l]<-gR[j]
	        	 l=l+1
	         }
        	 gResult1<-rbind(gResult1,gResult2)
         }
		 l<-1
         gResult2<-matrix("0",1,50)
         for(i in k[2*(length(k)/2-1)+1]:length(gR)){
         	 gResult2[1,l]<-gR[i]
         	 l<-l+1
         }
         gResult1<-rbind(gResult1,gResult2)
	 }else{
         l<-1
         gResult2<-matrix("0",1,50)
         for(i in k[2*(length(k)/2-1)+1]:length(gR)){
         	 gResult2[1,l]<-gR[i]
        	 l<-l+1
         }
         gResult1<-rbind(gResult1,gResult2)
     }
	 
		 gResult<-gResult1
		 gResult[,1]=as.character(gResult[,1])
		 list1<-list()		   
		 for(i in 1:dim(gResult)[1]){
		     k<-c()
		     k<-which(gResult[i,1]==gResult[,1])
		     list1<-c(list1,list(k=k))
		 }
		 data<-unique(list1)   
		 for(i in 1:length(data)){
		     if(length(data[[i]])==1){
			     gResult[data[[i]],1]=paste0(gResult[data[[i]],1],"_","1") 
			 }else{
			     for(j in 1:length(data[[i]])){
				     gResult[data[[i]][j],1]=paste0(gResult[data[[i]][j],1],"_",(data[[i]][j]-data[[i]][1]+1))
				 }
			 }
		 }
		 gResult2<-gResult
		 m<-c()
		 k<-0
         for(i in 10:dim(gResult2)[2]){
		     if(length(which(gResult2[,i]==0))==dim(gResult2)[1]){
			     m<-c(m,i)
			 }
		 }
         for(i in 1:length(m)){
		     gResult2<-gResult2[,-(m[i]-k)]
			 k=k+1
		 }		 
		 write.table(gResult2,"txtResult.txt",sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,quote=TRUE)
	     
	 if(exists("gResult2")&&dim(gResult2)[1]>0){
	     return(gResult2)
	 }else{
	     print("Cannot find subpathway!")
	 }
   }
}
