#' @title Download the network from GeneMania. 
#' @description
#'      SpidermiRdownload_net function will download the data 
#'
#' @param data The SpidermiRquery_spec_networks output
#' @examples
#' org<-SpidermiRquery_species(species)
#' net_shar_prot<-SpidermiRquery_spec_networks(organismID = org[9,],
#' network = "SHpd")
#' out_net<-SpidermiRdownload_net(data=net_shar_prot)
#' @export
#' @return Download GeneMania network 
SpidermiRdownload_net <- function(data){
    list_d<-list()
    for(i in 1:length(data)){
      print(paste("Downloading: ",data[i], " ... reference n. ", i, " of ", length(data), sep = ""))
      data<-unlist(data)
      list_d[[i]]<-read.table(data[i],header = TRUE,stringsAsFactors=FALSE)

    }
    return (list_d)  
}

#' @title Download human miRNA predicted database
#' @description SpidermiRdownload_miRNAprediction will download miRNA predicted target
#' @param mirna_list miRNA list of interest
#' @examples
#' mirna<-c('hsa-miR-567')
#' list<-SpidermiRdownload_miRNAprediction(mirna_list=mirna)
#' @export
#' @import stats
#' @import org.Hs.eg.db
#' @import miRNAtap.db
#' @importFrom miRNAtap getPredictedTargets 
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL2EG
#' @importFrom AnnotationDbi mappedkeys as.list
#' @return a dataframe with miRNA target validated interactions
SpidermiRdownload_miRNAprediction<-function(mirna_list){
  dop=list()
  dop2=list()
  for (k in  1:length(mirna_list)){
    print(paste("Processing...",mirna_list[k]))
    targets <- getPredictedTargets(mirna_list[k],species='hsa', method ='geom')
    if(is.null(targets)){
      dop2[[k]]<-mirna_list[k]
    }
    if(length(targets)!=0){
      dop[[k]]<-targets
      names(dop)[[k]]<-mirna_list[k]}
  }
  x <- org.Hs.egSYMBOL2EG
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  top <- matrix(0, length(xx), length(dop))
  rownames(top) <- names(xx)
  colnames(top)<- names(dop)
  #j=39
  #k=1
  for (j in  1:length(xx)){
    for (k in  1:length(dop)){
      if (length(intersect(xx[[j]],rownames(dop[[k]]))!=0)){
        #print (j)
        top[j,k]<-names(xx[j])
      }
    }  
  }
  top[top == 0] <- NA
  top[apply(top,1,function(x)any(!is.na(x))),]
  top<-t(top)
  top<-int(top)
  return(top)
}














#SpidermiRdownload_miRNAvalidate<-function(validated){

 #   site_mirtarbase<-.url_cache$get("miRTarBase")
  #  site_mirtarbase<-sub("s", "", site_mirtarbase)
   # test <- read.xls(site_mirtarbase, quote="",stringsAsFactors=FALSE)
    #test2<-test[test$X.Species..miRNA..=="\"Homo sapiens\"",]
    
    #pro<-as.data.frame(cbind(test2$X.miRNA.,test2$X.Target.Gene.))
    #dem<- as.data.frame(sapply(pro, function(x) gsub("\"", "", x)))

    
    #return(dem)
#}



#' @title Download miRNA validated database
#' @description SpidermiRdownload_miRNAextra_cir will download miRNA validated target
#' @param miRNAextra_cir parameter
#' @examples
#' list<-SpidermiRdownload_miRNAextra_cir(miRNAextra_cir)
#' @export
#' @import stats
#' @return a dataframe with miRNA target validated interactions
SpidermiRdownload_miRNAextra_cir<-function(miRNAextra_cir){
  # querying miRandola database (Extracellular Circulating microRNAs)
  site<-.url_cache$get("mirandola")
  mirandola<-read.delim(site,header = TRUE,quote = "",stringsAsFactors=FALSE)
  return(mirandola)
}





#' @title Download drug-gene interactions in  DGIdb 
#' @description SpidermiRdownload_drug_gene will download drug gene interactions
#' @param drug_gene parameter
#' @examples
#' drug_genetarget<-SpidermiRdownload_drug_gene(drug_gene)
#' @export
#' @import stats
#' @return a dataframe with miRNA target validated interactions
SpidermiRdownload_drug_gene<-function(drug_gene){
  # querying DGIdb
  site<-.url_cache$get("dgidb")
  dgidb<-read.delim(site,header = TRUE,quote = "",stringsAsFactors=FALSE)
  dgib_m<-dgidb[,c(1,8)]
  #url1<-.url_cache$get("matador")
  #con <- gzcon(url(url1))
  #txt <- readLines(con)
  #zz<-read.delim(textConnection(txt))
  #zz_m<-zz[,c(5,2)]
  #colnames(zz_m)<-colnames(dgib_m)
  dgib_m$gene_name<-toupper(dgib_m$gene_name)
  dgib_m$drug_name<-toupper(dgib_m$drug_name)
  #zz_m$gene_name<-toupper(zz_m$gene_name)
  #zz_m$drug_name<-toupper(zz_m$drug_name)
  #dgidb_matador<-rbind(zz_m,dgib_m)
  dgidb_matador<-dgib_m
  dgidb_matador_tot<-dgidb_matador[!duplicated(dgidb_matador), ]
  return(dgidb_matador_tot)
}




