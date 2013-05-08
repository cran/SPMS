getData <-
function(paths){
     library(XML)
     library(graph)
     library(RBGL)
	 setwd(paths)
	 if(length(which("nodes"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/nodes")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("nodes")
	 }
	 setwd(paths)
	 if(length(which("coords1"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/coords1")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("coords1")
	 }
	 setwd(paths)
	 if(length(which("depth1"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/depth1")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("depth1")
	 }
     setwd(paths)
	 if(length(which("broad1"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/broad1")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("broad1")
	 }
	 setwd(paths)
	 if(length(which("geneAndSymbol"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/geneAndSymbol")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("geneAndSymbol")
	 }
	 setwd(paths)
	 if(length(which("depth"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/depth")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("depth")
	 }
	 setwd(paths)
	 if(length(which("broad"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/broad")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("broad")
	 }
	 setwd(paths)
	 if(length(which("networkData"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/networkData")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("networkData")
	 }
     setwd(paths)
	 if(length(which("networkData1"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/networkData1")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("networkData1")
	 }
	 setwd(paths)
	 if(length(which("networkData2"==dir()))>0){
		 n<-sub(" ","",paste(paths,"/networkData2")) 
	     setwd(n)
		 unlink(dir())
	 }else{
		 dir.create("networkData2")
	 }
	 
	 n<-sub(" ","",paste(paths,"/xml")) 
	 setwd(n)
     for(p in 1:length(dir())){	     
		 n<-sub(" ","",paste(paths,"/xml")) 
	     setwd(n)
		 print(dir()[p])
		 pathway<-dir()[p]
		 name<-sub("xml","txt",pathway)
		 if(file.info(pathway)$size==2)next	 
		 
		 top<-xmlRoot(xmlTreeParse(pathway))		 
		 KOList1<-character()
         KOListIndex1<-0
	     KOList2<-character()
         KOListIndex2<-0
	     KOList_g<-character()
         for(i in 1:length(top)){
            if(xmlName(top[[i]])=="entry"){
                  if(!is.na(xmlAttrs(top[[i]])["type"])){
                        if(!is.na(xmlAttrs(top[[i]])["name"])){
                              if(xmlAttrs(top[[i]])["type"]=="gene"){
                                  tmp1<-unlist(xmlAttrs(top[[i]])["name"])
							      tmp2<-unlist(strsplit(xmlAttrs(top[[i]])["id"]," "))
                                  if(length(tmp1)>0){
                                     for(k in 1:length(tmp1)){
                                         KOListIndex1<-KOListIndex1+1
                                         KOList1[KOListIndex1]<-tmp1[k]
                                     }
                                  } 
							      if(length(tmp2)>0){
                                     for(k in 1:length(tmp2)){
                                         KOListIndex2<-KOListIndex2+1
                                         KOList2[KOListIndex2]<-tmp2[k]
                                     }
                                  }
                              }
                        }                       
                  }
            }
         }
	     for(i in 1:length(KOList1)){
	         KOList_g<-rbind(KOList_g,cbind(KOList2[i],KOList1[i]))
	     }
         gene_node<-KOList_g
	  
	     KOList3<-character()
         KOListIndex3<-0
	     KOList4<-character()
         KOListIndex4<-0
	     KOList_r<-character()
         for(i in 1:length(top)){
            if(xmlName(top[[i]])=="relation"){
                  if(!is.na(xmlAttrs(top[[i]])["type"])){
                        if(!is.na(xmlAttrs(top[[i]])["entry1"])){
                              if(xmlAttrs(top[[i]])["type"]=="ECrel"){
                                  tmp3<-unlist(strsplit(xmlAttrs(top[[i]])["entry1"]," "))
							      tmp4<-unlist(strsplit(xmlAttrs(top[[i]])["entry2"]," "))
                                  if(length(tmp3)>0){
                                      for(k in 1:length(tmp3)){
                                         KOListIndex3<-KOListIndex3+1
                                         KOList3[KOListIndex3]<-tmp3[k]
                                      }
                                 }
							      if(length(tmp4)>0){
                                      for(k in 1:length(tmp4)){
                                          KOListIndex4<-KOListIndex4+1
                                          KOList4[KOListIndex4]<-tmp4[k]
                                      }
                                  }
                              }
                        }                       
                  }
            }
         }
	  
	     for(i in 1:length(KOList3)){
	         KOList_r<-rbind(KOList_r,cbind(KOList3[i],KOList4[i]))		  
	     }	    
         node_relation<-KOList_r		 
		 
	  if(is.na(node_relation[1,1])==FALSE){
	  
	     
	     node<-c()
		 node<-unique(c(node_relation[,1],node_relation[,2]))
		 temp<-c()
		 for(i in 1:length(node)){
		     if(length(which(gene_node[,1]==node[i]))==0){
			     temp<-c(temp,node[i])
			 }			 
		 }
		 for(i in 1:length(temp)){
		     m<-which(node_relation[,1]==temp[i])
			 n<-which(node_relation[,2]==temp[i])
			 if(length(m)>0){
			     for(j in 1:length(m)){
				     node_relation[m[j],]<-"0"
				 }
			 }
			 if(length(n)>0){
			     for(j in 1:length(n)){
				     node_relation[n[j],]<-"0"
				 }
			 }
		 }
		 node_relation<-unique(node_relation)
		 for(i in 1:dim(node_relation)[1]){
		     if(node_relation[i,1]=="0"){
			     node_relation<-node_relation[-i,]
				 break
			 }			 
		 }
         if(dim(node_relation)[1]==0)next	     
         
		 k1<-sub(" ","",paste(paths,"/networkData"))
	     setwd(k1)
		 write.table(node_relation,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
		 
         node_re<-node_relation
		 for(i in 1:dim(node_re)[1]){
		     k<-0
		     k=which(node_re[i,1]==gene_node[,1])
			 node_re[i,1]=gene_node[k,2]
			 l<-0
		     l=which(node_re[i,2]==gene_node[,1])
			 node_re[i,2]=gene_node[l,2]
		 }
		 k2<-sub(" ","",paste(paths,"/networkData2"))
	     setwd(k2)
		 write.table(node_re,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
		 
         		 
	     KOList5<-character()
         KOListIndex5<-0
	     gene_symbol<-character()
         for(i in 1:length(top)){
            if(xmlName(top[[i]])=="entry"){
                  if(!is.na(xmlAttrs(top[[i]])["type"])){
                        if(!is.na(xmlAttrs(top[[i]])["name"])){
                              if(xmlAttrs(top[[i]])["type"]=="gene"){
							     tmp5<-unlist(xmlChildren(top[[i]]))							 
							     if(length(tmp5)>0){
                                     for(k in 1:length(tmp5)){
                                         KOListIndex5<-KOListIndex5+1
                                         KOList5[KOListIndex5]<-tmp5[k]
                                     }
                                 }
                              }
                        }                       
                  }
            }
         }
	     gene_symbol<-KOList5
         gene_AllSymbol<-gene_symbol
		 
	     gene_symbol<-c()
	     for(i in 1:length(gene_AllSymbol)){
		     if(gene_AllSymbol[i]=="graphics"){
			     gene_symbol<-c(gene_symbol,gene_AllSymbol[i+1])
			 }
	     }
		 
	     geneAndSymbol<-c()
	     nodes<-gene_node[,2]
         for(i in 1:dim(gene_node)[1]){
	         geneAndSymbol<-rbind(geneAndSymbol,cbind(nodes[i],gene_symbol[i]))
	     }
	     k<-sub(" ","",paste(paths,"/geneAndSymbol"))
	     setwd(k)
		 write.table(geneAndSymbol,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
	    
         node_re<-node_relation
		 for(i in 1:dim(node_re)[1]){
		     k<-0
		     k=which(node_re[i,1]==gene_node[,1])
			 m<-0
			 m=which(gene_node[k,2]==geneAndSymbol[,1])[1]
			 node_re[i,1]=geneAndSymbol[m,2]
			 l<-0
		     l=which(node_re[i,2]==gene_node[,1])
             n<-0
			 n=which(gene_node[l,2]==geneAndSymbol[,1])[1]
			 node_re[i,2]=geneAndSymbol[n,2]
		 }
		 k2<-sub(" ","",paste(paths,"/networkData1"))
	     setwd(k2)
		 write.table(node_re,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
		 
		
		 coord<-c()
		 for(i in 1:length(gene_AllSymbol)){
	         if(gene_AllSymbol[i]=="rectangle"){
			     gene_coord<-c()
				 gene_coord<-cbind(gene_coord,gene_AllSymbol[i+1],gene_AllSymbol[i+2])
				 coord<-rbind(coord,gene_coord)
			 }
	     }
		 
		 coords<-c()
		 coords_nodes<-gene_node[,1]
         for(i in 1:dim(coord)[1]){
	         coords<-rbind(coords,cbind(coords_nodes[i],coord[i,1],coord[i,2]))
	     }
	     k<-sub(" ","",paste(paths,"/coords1"))
	     setwd(k)
		 write.table(coords,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
		 
		 node_relation<-as.matrix(node_relation)
         nodes<-unique(c(node_relation[,1],node_relation[,2]))
         g1<-new("graphNEL",nodes=nodes,edgemode="undirected")     
         for(i in 1:dim(node_relation)[1]){
		     if(!isAdjacent(g1,node_relation[i,1],node_relation[i,2])){				
                 g2<-addEdge(node_relation[i,1],node_relation[i,2],g1,1)
                 g1<-g2
             }					
         }
         keggPathway<-g1
	 
         #node_relation<-as.matrix(node_relation)
         nodes<-unique(c(node_relation[,1],node_relation[,2]))
		 pathway_node<-c()
		 for(i in 1:length(nodes)){
		     k<-0
			 l<-0
			 k<-which(gene_node[,1]==nodes[i])
			 pathway_node<-c(pathway_node,geneAndSymbol[k,2])			 
		 }
	     k<-sub(" ","",paste(paths,"/nodes"))
	     setwd(k)
		 write.table(pathway_node,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
	     
		 
         edges(keggPathway)	
		 
         if(length(connComp(keggPathway))==1){
             col2<-list()		 
		     for(i in 1:length(nodes)){ 
			     col1<-c()
                 col1<-c(col1,(dfs(keggPathway,nodes[i],checkConn=TRUE))[[1]])
				 col2<-c(col2,list(col1=col1))
                 shortest_dfs<-c()
				 shortest_dfs<-c(shortest_dfs,0)
				 for(l in 2:length(col1)){
				     m<-c()						 
					 m<-sp.between(keggPathway,col1[1],col1[l], detail=FALSE)
					 shortest_dfs<-c(shortest_dfs,m[[1]]$length)						 
				 }
			     col2<-c(col2,list(shortest_dfs=shortest_dfs))				 
             }			 
		     depth<-col2
			 depth2<-matrix("0",length(nodes),length(depth))
			 for(l in 1:length(depth)){
			     if(l%%2==1){
			         for(o in 1:length(depth[[l]])){
			             depth2[o,l]<-as.character(depth[[l]][o])
				     }
				 }else{
				     for(o in 1:length(depth[[l]])){
			             depth2[o,l]<-as.character(depth[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/depth1"))
	         setwd(p)
		     write.table(depth2,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
	    	
	         for(i in 1:(length(depth)/2)){
                 for(j in 1:length(depth[[(2*i-1)]])){			 
	                     
						 k=which(depth[[(2*i-1)]][j]==gene_node[,1])
		                 depth[[(2*i-1)]][j]=gene_node[k,2]		  
		                 l=which(gene_node[k,2]	==geneAndSymbol[,1])[1]
						 depth[[(2*i-1)]][j]=geneAndSymbol[l,2]
                 }					 
             } 
             depth1<-matrix("0",length(nodes),length(depth))
			 for(l in 1:length(depth)){
			     if(l%%2==1){
			         for(o in 1:length(depth[[l]])){
			             depth1[o,l]<-as.character(depth[[l]][o])
				     }
				 }else{
				     for(o in 1:length(depth[[l]])){
			             depth1[o,l]<-as.character(depth[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/depth"))
	         setwd(p)
		     write.table(depth1,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
	    			 
	         list1<-list()
			 col<-list()
             for(i in 1:length(nodes)){
		         col1<-c()
                 col1<-c(col1,bfs(keggPathway,nodes[i],checkConn=TRUE))
                 list1<-c(list1,list(col1=col1))
                 shortest_bfs<-c()				 
				 shortest_bfs<-c(shortest_bfs,0)
				 for(l in 2:length(col1)){
				     m<-c()						 
					 m<-sp.between(keggPathway,col1[1],col1[l], detail=FALSE)
					 shortest_bfs<-c(shortest_bfs,m[[1]]$length)						 
				 }
			      list1<-c(list1,list(shortest_bfs=shortest_bfs))	
             }	
             broad<-list1
			 broad2<-matrix("0",length(nodes),length(broad))
			 for(l in 1:length(broad)){
			     if(l%%2==1){
			         for(o in 1:length(broad[[l]])){
			             broad2[o,l]<-as.character(broad[[l]][o])
				     }
				 }else{
				     for(o in 1:length(broad[[l]])){
			             broad2[o,l]<-as.character(broad[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/broad1"))
	         setwd(p)
		     write.table(broad2,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)

			 for(i in 1:(length(broad)/2)){
                 for(j in 1:length(broad[[(2*i-1)]])){			 
	                     
						 k=which(broad[[(2*i-1)]][j]==gene_node[,1])
		                 broad[[(2*i-1)]][j]=gene_node[k,2]		  
		                 l=which(gene_node[k,2]	==geneAndSymbol[,1])[1]
						 broad[[(2*i-1)]][j]=geneAndSymbol[l,2]
                 }					 
             } 
             broad1<-matrix("0",length(nodes),length(broad))
			 for(l in 1:length(broad)){
			     if(l%%2==1){
			         for(o in 1:length(broad[[l]])){
			             broad1[o,l]<-as.character(broad[[l]][o])
				     }
				 }else{
				     for(o in 1:length(broad[[l]])){
			             broad1[o,l]<-as.character(broad[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/broad"))
	         setwd(p)
		     write.table(broad1,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE)
		 }else{
		 	 col3<-c()
			 n<-0
			 depth<-list()
			 broad<-list()
			 col<-list()
			 list2<-list()
		     k<-length(connComp(keggPathway))
			 for(i in 1:k){
			     nodes_n<-connComp(keggPathway)[[i]]
				 col3<-c(col3,length(nodes_n))
				 n<-n+col3[i]
				 node_relation1<-character()
				 for(p in 1:dim(node_relation)[1]){
				     if(is.na(pmatch(node_relation[p,1],nodes_n))==FALSE&&is.na(pmatch(node_relation[p,2],nodes_n))==FALSE){
					     node_relation1<-rbind(node_relation1,cbind(as.character(node_relation[p,1]),as.character(node_relation[p,2])))
					 }
				 }
			     g1<-new("graphNEL",nodes=nodes_n,edgemode="undirected")     
                 for(j in 1:dim(node_relation1)[1]){
		             if(!isAdjacent(g1,node_relation1[j,1],node_relation1[j,2])){				
                         g2<-addEdge(node_relation1[j,1],node_relation1[j,2],g1,1)
                         g1<-g2
                     }					
                 }
                 keggPathway1<-g1
			     for(j in 1:length(nodes_n)){
				     col1<-c()
					 shortest_dfs<-c()
			         col1<-c(col1,(dfs(keggPathway1,nodes_n[j],checkConn=TRUE))[[1]])
					 col<-c(col,list(col1=col1))
					 shortest_dfs<-c(shortest_dfs,0)
					 for(l in 2:length(col1)){
					     m<-c()						 
						 m<-sp.between(keggPathway1,col1[1],col1[l], detail=FALSE)
						 shortest_dfs<-c(shortest_dfs,m)						 
					 }
					 col<-c(col,list(shortest_dfs=shortest_dfs))
					 col2<-c()
					 shortest_bfs<-c()
			         col2<-c(col2,bfs(keggPathway1,nodes_n[j],checkConn=TRUE))
					 list2<-c(list2,list(col2=col2))
					 shortest_bfs<-c(shortest_bfs,0)
					 for(l in 2:length(col2)){
					     m<-c()						 
						 m<-sp.between(keggPathway1,col2[1],col2[l], detail=FALSE)
						 shortest_bfs<-c(shortest_bfs,m)						 
					 }
					 list2<-c(list2,list(shortest_bfs=shortest_bfs))
			     }                 				 
             }
             m<-max(col3)             			 
             depth<-col
			 depth2<-matrix("0",m,n*2)
			 for(l in 1:length(depth)){
			     if(l%%2==1){
			         for(o in 1:length(depth[[l]])){
			             depth2[o,l]<-as.character(depth[[l]][o])
				     }
				 }else{
				     for(o in 1:length(depth[[l]])){
			             depth2[o,l]<-as.character(depth[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/depth1"))
	         setwd(p)
		     write.table(depth2,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,na="NA")	    
		 
			 for(i in 1:(length(depth)/2)){
                 for(j in 1:length(depth[[(2*i-1)]])){			 	                     
						 k=which(depth[[(2*i-1)]][j]==gene_node[,1])
		                 if(length(k)>0){
					         depth[[(2*i-1)]][j]=gene_node[k,2]		  
		                     l=which(gene_node[k,2]	==geneAndSymbol[,1])[1]
						     depth[[(2*i-1)]][j]=geneAndSymbol[l,2]
						 }
                 }					 
             }             			 
			 depth1<-matrix("0",m,n*2)
			 for(l in 1:length(depth)){
			     if(l%%2==1){
			         for(o in 1:length(depth[[l]])){
			             depth1[o,l]<-as.character(depth[[l]][o])
				     }
				 }else{
				     for(o in 1:length(depth[[l]])){
			             depth1[o,l]<-as.character(depth[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/depth"))
	         setwd(p)
		     write.table(depth1,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,na="NA")	    
		 			 
             broad<-list2
			 broad2<-matrix("0",m,n*2)
			 for(l in 1:length(broad)){
			     if(l%%2==1){
			         for(o in 1:length(broad[[l]])){
			             broad2[o,l]<-as.character(broad[[l]][o])
				     }
				 }else{
				     for(o in 1:length(broad[[l]])){
			             broad2[o,l]<-as.character(broad[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/broad1"))
	         setwd(p)
		     write.table(broad2,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,na="NA")	    
			 
             for(i in 1:(length(broad)/2)){
                 for(j in 1:length(broad[[(2*i-1)]])){			 	                     
						 k=which(broad[[(2*i-1)]][j]==gene_node[,1])
		                 if(length(k)>0){
					         broad[[(2*i-1)]][j]=gene_node[k,2]		  
		                     l=which(gene_node[k,2]	==geneAndSymbol[,1])[1]
						     broad[[(2*i-1)]][j]=geneAndSymbol[l,2]
						 }
                 }					 
             }             			 
			 broad1<-matrix("0",m,n*2)
			 for(l in 1:length(broad)){
			     if(l%%2==1){
			         for(o in 1:length(broad[[l]])){
			             broad1[o,l]<-as.character(broad[[l]][o])
				     }
				 }else{
				     for(o in 1:length(broad[[l]])){
			             broad1[o,l]<-as.character(broad[[l]][o][[1]])
				     }
				 }
			 }
	         p<-sub(" ","",paste(paths,"/broad"))
	         setwd(p)
		     write.table(broad1,name,quote=T,sep="\t",col.names=FALSE,append=TRUE,row.names=FALSE,na="NA")	    
		  }
	   }
    } 
	 return()
}
