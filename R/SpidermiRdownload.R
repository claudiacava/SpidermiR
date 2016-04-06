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
      data<-unlist(data)
      list_d[[i]]<-read.table(data[i],header = TRUE,stringsAsFactors=FALSE)
      #print(paste(data[i], " ... reference n. ", i, " of ", length(data), sep = ""))
    }
    return (list_d)  
}

#' @title Download both miRNA target and the gene-drug interaction from PharmacomiR database
#' @description SpidermiRdownload_pharmacomir will download miRNA Pharmacogenomic data
#' @param pharmacomir variable
#' @export
#' @return a dataframe with gene-drug, and miR-gene associations
#' @examples
#' mir_pharmaco<-SpidermiRdownload_pharmacomir(pharmacomir)
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

