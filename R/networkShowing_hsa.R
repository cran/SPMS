networkShowing_hsa <-
function(paths,pathway){
     library(igraph)
     n<-sub(" ","",paste(paths,"/networkData2")) 
	 setwd(n)
	 name<-gsub(" ","",paste(pathway,".txt"))
	 data<-read.table(name,sep="\t",header=FALSE)
	 d=data.frame(p1=data[,1],p2=data[,2])
	 g = graph.data.frame(d, directed = FALSE)
	 tkplot(g,vertex.label=V(g)$name,vertex.size=20)
}
