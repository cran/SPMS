ErgodicPathShowing <-
function(paths,dORb,pathway,q){
     library(igraph)
	 if(dORb=="depth"){
	     n<-sub(" ","",paste(paths,"/depth"))
	 }else{
	     n<-sub(" ","",paste(paths,"/broad"))
	 }
	 setwd(n)
	 name<-gsub(" ","",paste(pathway,".txt"))
	 entry_id_depth<-read.table(name,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 n1<-sub(" ","",paste(paths,"/temp/resulttemp"))
	 setwd(n1)
	 name1<-gsub(" ","",paste(pathway,".txt"))
	 resulttemp<-read.table(name1,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 
	 ergodic<-entry_id_depth[,resulttemp[q,1]]
	 data<-c()
	 for(i in 1:(length(ergodic)-1)){
	     col<-c()
	     col<-cbind(ergodic[i],ergodic[i+1])
		 data<-rbind(data,col)
	 }
	 
	 n2<-sub(" ","",paste(paths,"/networkData1")) 
	 setwd(n2)
	 name<-gsub(" ","",paste(pathway,".txt"))
	 data1<-read.table(name,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 data2<-matrix(0.1,dim(data1)[1],3)
	 data2[,1]<-data1[,1]
	 data2[,2]<-data1[,2]
	 data1<-data2
	 
	 for(i in 1:dim(data)[1]){
	     k<-c()
		 k=which(data[i,1]==data2[,1])
		 if(length(k)>0){
		     for(j in 1:length(k)){
		         if(data[i,2]==data2[k[j],2]){
	                 data2[k[j],3]=5		 
		    	 }else{
				     data2<-rbind(data2,cbind(data[i,1],data[i,2],2.5))
				 }
		     }
		 }else{
			 data2<-rbind(data2,cbind(data[i,1],data[i,2],2.5))
		 }
		 l<-c()
		 l=which(data[i,2]==data2[,1])
		 if(length(l)>0){
		     for(j in 1:length(l)){
		         if(data[i,1]==data2[l[j],2]){
	                 data2[l[j],3]=5		 
		    	 }else{
				     data2<-rbind(data2,cbind(data[i,1],data[i,2],2.5))
				 }
		     }
		 }else{
			 data2<-rbind(data2,cbind(data[i,1],data[i,2],2.5))
		 }
	 }
	 data3<-rbind(data2[1:dim(data1)[1],],unique(data2[dim(data1)[1]:dim(data2)[1],]))

	 for(i in 1:dim(data1)[1]){
	     for(j in (dim(data1)[1]+1):dim(data3)[1]){
		     if((data3[i,1]==data3[j,1]&&data3[i,2]==data3[j,2])||(data3[i,1]==data3[j,2]&&data3[i,2]==data3[j,1])){
			     data3[j,]=0
			 }
		 }
	 }
	 data3<-unique(data3)
	 for(i in 1:dim(data3)[1]){
	     if(data3[i,1]==0){
		     data3<-data3[-i,]
			 break
		 }
	 }
	 
	 d=data.frame(p1=data3[,1],p2=data3[,2])
	 g = graph.data.frame(d, directed = FALSE)
	 tkplot(g,vertex.size=20,edge.width=as.numeric(data3[,3]))
	 
}
