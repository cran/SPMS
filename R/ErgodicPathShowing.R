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
	 gene_ergodic_data<-read.table(name,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 if(dORb=="depth"){
	     n<-sub(" ","",paste(paths,"/depth1"))
	 }else{
	     n<-sub(" ","",paste(paths,"/broad1"))
	 }
	 setwd(n)
	 name<-gsub(" ","",paste(pathway,".txt"))
	 id_ergodic_data<-read.table(name,sep="\t",header=FALSE,stringsAsFactors=FALSE)
#====================input entry_id network ergodic_data=====================	  
	 n<-sub(" ","",paste(paths,"/networkData")) 
	 setwd(n)
	 id_net<-read.table(name,sep="\t",header=FALSE,stringsAsFactors=FALSE)
#=====================input entry_id ergodic ergodic_data=======================

	 n1<-sub(" ","",paste(paths,"/temp/resulttemp"))
	 setwd(n1)
	 name1<-gsub(" ","",paste(pathway,".txt"))
	 resulttemp<-read.table(name1,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 
	 id_ergodic<-id_ergodic_data[,resulttemp[q,1]]
	 ergodic_id_data<-c()
	 for(i in 1:(length(id_ergodic)-1)){
	     col<-c()
	     col<-cbind(id_ergodic[i],id_ergodic[i+1])
		 ergodic_id_data<-rbind(ergodic_id_data,col)
	 }
	 	 
	 ergodic<-gene_ergodic_data[,resulttemp[q,1]]
	 ergodic_data<-c()
	 for(i in 1:(length(ergodic)-1)){
	     col<-c()
	     col<-cbind(ergodic[i],ergodic[i+1])
		 ergodic_data<-rbind(ergodic_data,col)
	 }
	 
#=======================input netwotk data=====================	 
	 n2<-sub(" ","",paste(paths,"/networkData1")) 
	 setwd(n2)
	 name<-gsub(" ","",paste(pathway,".txt"))
	 gene_net<-read.table(name,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 gene_data<-matrix(0.1,(dim(gene_net)[1]+dim(ergodic_data)[1]),3)
	 gene_data[,1]<-c(gene_net[,1],ergodic_data[,1])
	 gene_data[,2]<-c(gene_net[,2],ergodic_data[,2])
	 gene_data[(dim(gene_net)[1]+1):dim(gene_data)[1],3]=3
	 
	 id_data<-matrix(0.1,(dim(id_net)[1]+dim(ergodic_data)[1]),2)
	 id_data[,1]<-c(id_net[,1],ergodic_id_data[,1])
	 id_data[,2]<-c(id_net[,2],ergodic_id_data[,2])
	 
     data<-c(gene_data[,1],gene_data[,2])
	 data1<-c(id_data[,1],id_data[,2])
	 
	 for(i in 1:length(data)){
		 k=which(data[i]==data)
		 if(length(k)>1){
		     for(j in 1:(length(k)-1)){
		         for(l in (j+1):length(k)){
		             if(data1[k[j]]!=data1[k[l]]){
		    		     data[k[l]]=paste(data[k[l]],"")					 
	    			 }			     
	    		 }
				 break
		     }
		 }		 
	 }
	 data3<-cbind(data[1:(length(data)/2)],data[(length(data)/2+1):length(data)],gene_data[,3])
	 
	 for(i in 1:dim(gene_net)[1]){
	     for(j in (dim(gene_net)[1]+1):dim(data3)[1]){
		     if(data3[i,1]==data3[j,1]&&data3[i,2]==data3[j,2]||data3[i,1]==data3[j,2]&&data3[i,2]==data3[j,1]){
			     data3[j,]=0
				 data3[i,3]=3
			 }
		 }
	 }
	 data3<-unique(data3)
	 k=which("0"==data3[,2])
	 if(length(k)>0){
	     data3<-data3[-k,]
	 }
	 
	 d=data.frame(p1=data3[,1],p2=data3[,2])
	 g = graph.data.frame(d, directed = FALSE)
	 tkplot(g,vertex.size=20,edge.width=as.numeric(data3[,3]))
	 
}
