#' @title Download the network from GeneMania. 
#' @description
#'      SpidermiRdownload_net function will download the data 
#'
#' @param data The SpidermiRquery_spec_networks output
#' @examples
#' org<-SpidermiRquery_species(species)
#' net_shar_prot<-SpidermiRquery_spec_networks(organismID = org[9,],
#'                                        network = "SHpd")
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

#' @title Download both miRNA target and the gene-drug interaction from PharmacomiR database
#' @description SpidermiRdownload_pharmacomir will download miRNA Pharmacogenomic data
#' @param pharmacomir variable
#' @examples
#' mir_pharmaco<-SpidermiRdownload_pharmacomir(pharmacomir=pharmacomir)
#' @export
#' @import stats
#' @return a dataframe with gene-drug, and miR-gene associations
SpidermiRdownload_pharmacomir<-function(pharmacomir){
  # querying Pharmaco-miR database (Pharmaco-miR validated interaction)
  pharm_miR <- .url_cache$get("pharmacomir")
  pharm_miR_db<-read.csv(pharm_miR,header = TRUE,stringsAsFactors=FALSE)
  pharm_miR_db$miRNA <- as.character(sub("miR-","hsa-miR-", pharm_miR_db$miRNA))
  pharm_miR_db$miRNA <- as.character(sub("let-","hsa-let-", pharm_miR_db$miRNA))
  pharm_miR_db$Drug<-as.character(sub("3,3'-","", pharm_miR_db$Drug))
  pharm_miR_db$Drug<-as.character(sub("5-","", pharm_miR_db$Drug))
  pmir<-cbind(pharm_miR_db$miRNA,pharm_miR_db$Gene)
  pmir<-as.data.frame(pmir[!duplicated(pmir), ]) 
  pmir2<-cbind(pharm_miR_db$Gene,pharm_miR_db$Drug)
  pmir2<-as.data.frame(pmir2[!duplicated(pmir2), ]) 
  pmir_c<-rbind(pmir,pmir2)
  i <- sapply(pmir_c, is.factor)
  pmir_c[i] <- lapply(pmir_c[i], as.character)
  return(pmir_c)
  #SpidermiRvisualize_mirnanet(data=pmir_c)
}


#' @title Download miRNA predicted database
#' @description SpidermiRdownload_miRNAprediction will download miRNA predicted target
#' @param mirna_list miRNA list of interest
#' @examples
#' mirna<-c('hsa-miR-567','hsa-miR-566')
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



#' @title Download miRNA validated database
#' @description SpidermiRdownload_miRNAprediction will download miRNA validated target
#' @param validated parameter
#' @examples
#' list<-SpidermiRdownload_miRNAvalidate(validated)
#' @export
#' @import stats
#' @return a dataframe with miRNA target validated interactions
SpidermiRdownload_miRNAvalidate<-function(validated){
    # querying miRtar database (validated interaction miRNA-gene)
    site_mir2disease<-.url_cache$get("miRtar")
    mir2disease<-read.delim(site_mir2disease,header = FALSE,quote = "",stringsAsFactors=FALSE)
    # querying miRNA WALK database (validated interaction miRNA-gene)
    temp <- tempfile()
    download.file(.url_cache$get("miRwalk"),temp)
    sx<-unz(temp,"hsa-vtm-gene.rdata.Rdata")
    load(sx)
    id<-t(sapply(id, '[', 1:max(sapply(id, length)))) 
    se<-int(id)
    # merging miRtar and miRNA walk information
    mir2disease$V3<-NULL
    mir2disease$V4<-NULL
    mir_validated_targe<-rbind(mir2disease,se)
    return(mir_validated_targe)
}



#' @title Download miRNA validated database
#' @description SpidermiRdownload_miRNAprediction will download miRNA validated target
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




