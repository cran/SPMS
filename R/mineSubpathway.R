mineSubpathway <-
function(paths,interestOfGene,dORb,name,inG,p_value,broad_path){
	 fPathway<-fPathway(paths,name)
	 if(name=="hsa"||substr(name,1,3)=="hsa"){
	     num<-39114
	 }else{
	     num<-getNodesNum(paths)
	 }
	 # interestOfGene<-inputInterestGene()
	 if(nchar(name)<5){
	     if(dORb=="depth"){
             data<-depth_subPathwayMinner(paths,num,name,interestOfGene,p_value)
	         if(data!="Cannot find subpathway!"&&dim(data)[1]>1){
			     data1<-depth_zhenghe(paths,dORb,inG)	
			 }
         }else{
             data<-broad_subPathwayMinner(paths,num,name,interestOfGene,p_value,broad_path)
	         if(data!="Cannot find subpathway!"&&dim(data)[1]>1){
			     data1<-broad_zhenghe(paths,dORb,inG)	
			 }
		 }			 
	 }else{
	    if(dORb=="depth"){
             data<-depth_subPathwayMinner(paths,num,name,interestOfGene,p_value)
	         if(data!="Cannot find subpathway!"&&dim(data)[1]>1){
			     data1<-depth_zhenghe(paths,dORb,inG)	
			 }	
         }else{
             data<-broad_subPathwayMinner(paths,num,name,interestOfGene,p_value,broad_path)
	         if(data!="Cannot find subpathway!"&&data=="Cannot find subpathway!"){}
			 else{
			     data1<-broad_zhenghe(paths,dORb,inG)	
			 }
		 }	
	 }
}
