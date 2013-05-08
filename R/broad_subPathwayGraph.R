broad_subPathwayGraph <-
function(paths,q,pathway){
	 q<-as.numeric(q)
	 name<-pathway
     pathway<-sub(" ",".",paste(pathway,"txt"))
     m<-sub(" ","",paste(paths,"/geneAndSymbol"))
	 setwd(m)
	 geneAndSymbol<-read.table(pathway,sep="",header=FALSE,stringsAsFactors=FALSE)
	 n<-sub(" ","",paste(paths,"/temp/resulttemp"))
	 setwd(n)
	 gSymbol<-read.delim(pathway,header=FALSE,sep="\t",na.strings="NA",stringsAsFactors=FALSE)
	 n1<-sub(" ","",paste(paths,"/txtResult"))
	 setwd(n1)
	 gR<-read.table("txtResult.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
	 
     gResult<-c()
     for(i in 1:dim(gR)[1]){
    	 if(substr(gR[i,1],1,nchar(name))==name){
    		 gResult<-rbind(gResult,gR[i,])
         }
     }
	 
	 n2<-sub(" ","",paste(paths,"/broad1"))
	 setwd(n2)
	 broad1<-read.table(pathway,header=FALSE,sep="\t",na.strings="NA",stringsAsFactors=FALSE)
	 n3<-sub(" ","",paste(paths,"/broad"))
	 setwd(n3)
	 broad<-read.table(pathway,header=FALSE,sep="\t",na.strings="NA",stringsAsFactors=FALSE)
	 n4<-sub(" ","",paste(paths,"/coords"))
	 setwd(n4)
	 coords<-read.table(pathway,header=FALSE,sep="\t",na.strings="NA",stringsAsFactors=FALSE)

	 data(ImageSize,envir=environment())
	 ImageSize<-ImageSize
	 
     library(png)
	 library(EBImage)
	 o<-sub(" ","",paste(paths,"/picture"))
	 setwd(o)
     pdffile<-paste0(gResult[q,1],".pdf")
     pdf(file=pdffile, bg="transparent")

     op <- par(bg = "thistle")
     name1<-substr(name,(nchar(name)-4),nchar(name))
	 name1<-paste0("a",name1)
     a<-which(name1==ImageSize[,1])
     plot(c(0, as.numeric(as.character(ImageSize[a,2]))), c(0, as.numeric(as.character(ImageSize[a,3]))), type = "n", xlab="", ylab="",axes=FALSE,ann=FALSE)
	 p<-nchar(name)
	 species<-c()
	 if(p==8){
	     species<-substr(name,1,3)
	 }else if(p==9){
	     species<-substr(name,1,4)
	 }
	 # url=paste0("http://www.genome.jp/kegg/pathway/",species,"/",name,".png")
	 # image <-readImage(url)
	 p<-sub(" ","",paste(paths,"/png"))
	 setwd(p)
	 name1<-sub(" ",".",paste(name,"png"))
     image <- readPNG(name1,native=TRUE)
     rasterImage(image, 0, 0, as.numeric(as.character(ImageSize[a,2])), as.numeric(as.character(ImageSize[a,3])), interpolate=FALSE)
     par(new=TRUE)
	 for(i in gSymbol[q,2]:gSymbol[q,3]){
	     k<-0
		 k=which(coords[,1]==broad1[i,gSymbol[q,1]])		 
		 rect((as.numeric(coords[k,2])-23),abs((as.numeric(coords[k,3])+9)-as.numeric(as.character(ImageSize[a,3]))),(as.numeric(coords[k,2])+23),abs((as.numeric(coords[k,3])-9)-as.numeric(as.character(ImageSize[a,3]))), border="#FF00FF")
     }
	 for(i in 10:(9+gResult[q,4])){
	     for(j in gSymbol[q,2]:gSymbol[q,3]){
	         m<-strsplit(broad[j,gSymbol[q,1]],", ")[[1]]
    		 if(is.element(gResult[q,i],m)==TRUE){
			     k<-0
				 k=which(coords[,1]==broad1[j,gSymbol[q,1]])
    		     rect((as.numeric(coords[k,2])-23),abs((as.numeric(coords[k,3])+9)-as.numeric(as.character(ImageSize[a,3]))),(as.numeric(coords[k,2])+23),abs((as.numeric(coords[k,3])-9)-as.numeric(as.character(ImageSize[a,3]))), border="#FF0000")
             }
    	 }
	 }	
   dev.off()
   getwd()
}
