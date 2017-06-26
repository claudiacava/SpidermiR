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



#' @title Download miRNA validated database
#' @description SpidermiRdownload_miRNAprediction will download miRNA validated target
#' @param validated parameter
#' @examples
#' list<-SpidermiRdownload_miRNAvalidate(validated)
#' @export
#' @import stats
#' @importFrom gdata read.xls 
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
    site_mirtarbase<-.url_cache$get("miRTarBase")
    test <- read.xls(site_mirtarbase, quote="",stringsAsFactors=FALSE)
    pro<-as.data.frame(cbind(test$X.miRNA.,test$X.Target.Gene.))
    dem<- print(as.data.frame(sapply(pro, function(x) gsub("\"", "", x))))
    mir_validated_targe3<-rbind(mir_validated_targe,dem)
    return(mir_validated_targe3)
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





#' @title Download Long Non-Coding RNA (LNC-RNA) from  lncRNome database
#' @description SpidermiRdownload_LNC-RNA will download LNC-RNA 
#' @param miRNALNC_RNA parameter
#' @examples
#' list<-SpidermiRdownload_miRNALNC_RNA(miRNALNC_RNA)
#' @export
#' @import stats
#' @return a dataframe with miRNALNC_RNA
SpidermiRdownload_miRNALNC_RNA<-function(miRNALNC_RNA){
  # querying miRandola database (Extracellular Circulating microRNAs)
  site<-.url_cache$get("LNC")
  LNC<-read.delim(site,header = TRUE,quote = "",stringsAsFactors=FALSE)
  return(LNC)
}




#' @title Download microRNAs binding sites on long non coding RNA from  lncRNome database
#' @description SpidermiRdownload_miRNALNC_miRNA-RNA will download microRNAs binding sites on long non coding RNA
#' @param miRNALNC_miRNA parameter
#' @examples
#' list_LNC_miRNA<-SpidermiRdownload_miRNALNC_miRNA(miRNALNC_miRNA)
#' @export
#' @import stats
#' @return a dataframe with miRNALNC_miRNA
SpidermiRdownload_miRNALNC_miRNA<-function(miRNALNC_miRNA){
  site<-.url_cache$get("LNC_mirna")
  LNC_miRNA<-read.delim(site,header = TRUE,quote = "",stringsAsFactors=FALSE)
  return(LNC_miRNA)
}









#' @title Download miRNA tissue specific from Tiger Database
#' @description SpidermiRdownload_miR_tissue_spec will download miRNAs and genes tissue specific according to Tiger database (Liu X et al. 2008, Bioinformatics) 
#' @param miRNATS 
#' @examples
#' list_miRNA_TS<-SpidermiRdownload_miR_tissue_spec(miRNATS)
#' @export
#' @import stats
#' @return a dataframe with gene and miRNA tissue specific 
SpidermiRdownload_miR_tissue_spec<-function(miRNATS){
  # querying TIGER database
  site<-.url_cache$get("Tiger")
  nc<-max(count.fields(site, sep=""))
  tiger_db<-read.table(site,sep="",col.names=paste("v",1:nc,sep="."),fill=T,stringsAsFactors=FALSE)
  tiger_db<- matrix(unlist(tiger_db),ncol=nc)
  
  site_mapping<-.url_cache$get("Tiger_mapping")
  nc<-max(count.fields(site_mapping, sep=""))
  x<-read.table(site_mapping,sep="",col.names=paste("v",1:nc,sep="."),fill=T,stringsAsFactors=FALSE)
  tiger_mapping<-   matrix(unlist(x),ncol=nc)
  
  
  am<-match(tiger_db[,1],tiger_mapping[,2])
  st3 <- cbind( tiger_db, tiger_mapping[am,] )
  st3_final<-st3[complete.cases(st3),]
  colnames(st3_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am2<-match(tiger_db[,1],tiger_mapping[,3])
  st4 <- cbind( tiger_db, tiger_mapping[am2,] )
  st4_final<-st4[complete.cases(st4),]
  colnames(st4_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am3<-match(tiger_db[,1],tiger_mapping[,4])
  st5 <- cbind( tiger_db, tiger_mapping[am3,] )
  st5_final<-st5[complete.cases(st5),]
  colnames(st5_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am4<-match(tiger_db[,1],tiger_mapping[,5])
  st6 <- cbind( tiger_db, tiger_mapping[am4,] )
  st6_final<-st6[complete.cases(st6),]
  colnames(st6_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am5<-match(tiger_db[,1],tiger_mapping[,6])
  st7 <- cbind( tiger_db, tiger_mapping[am5,] )
  st7_final<-st7[complete.cases(st7),]
  colnames(st7_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am6<-match(tiger_db[,1],tiger_mapping[,7])
  st8 <- cbind( tiger_db, tiger_mapping[am6,] )
  st8_final<-as.matrix(t(st8[complete.cases(st8),]))
  #colnames(st8_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am7<-match(tiger_db[,1],tiger_mapping[,8])
  st9 <- cbind( tiger_db, tiger_mapping[am7,] )
  st9_final<-as.matrix(t(st9[complete.cases(st9),]))
  #colnames(st9_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  
  am8<-match(tiger_db[,1],tiger_mapping[,9])
  st10 <- cbind( tiger_db, tiger_mapping[am8,] )
  st10_final<-st10[complete.cases(st10),]
  colnames(st10_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  
  
  am9<-match(tiger_db[,1],tiger_mapping[,10])
  st11 <- cbind( tiger_db, tiger_mapping[am9,] )
  st11_final<-st11[complete.cases(st11),]
  colnames(st11_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am10<-match(tiger_db[,1],tiger_mapping[,11])
  st12 <- cbind( tiger_db, tiger_mapping[am10,] )
  st12_final<-st12[complete.cases(st12),]
  colnames(st12_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am11<-match(tiger_db[,1],tiger_mapping[,12])
  st13 <- cbind( tiger_db, tiger_mapping[am11,] )
  st13_final<-st13[complete.cases(st13),]
  colnames(st13_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  am12<-match(tiger_db[,1],tiger_mapping[,13])
  st14 <- cbind( tiger_db, tiger_mapping[am12,] )
  st14_final<-st14[complete.cases(st14),]
  colnames(st14_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  
  am13<-match(tiger_db[,1],tiger_mapping[,14])
  st15 <- cbind( tiger_db, tiger_mapping[am13,] )
  st15_final<-st15[complete.cases(st15),]
  colnames(st15_final) <- c("unigene1", "tissue1","tissue2","tissue3","tissue4","gene","unigene2","unigene3","unigene4","unigene5","unigene6","unigene7","unigene8","unigene9","unigene10","unigene11","unigene12","unigene13","unigene14")
  
  
  
  
  final_tiger<-cbind(st3_final[,6],st3_final[,2],st3_final[,3],st3_final[,4],st3_final[,5])
  final_tiger2<-cbind(st4_final[,6],st4_final[,2],st4_final[,3],st4_final[,4],st4_final[,5])
  final_tiger3<-cbind(st5_final[,6],st5_final[,2],st5_final[,3],st5_final[,4],st5_final[,5])
  final_tiger4<-cbind(st6_final[,6],st6_final[,2],st6_final[,3],st6_final[,4],st6_final[,5])
  final_tiger5<-cbind(st7_final[,6],st7_final[,2],st7_final[,3],st7_final[,4],st7_final[,5])
  final_tiger6<-cbind(st8_final[,6],st8_final[,2],st8_final[,3],st8_final[,4],st8_final[,5])
  final_tiger7<-cbind(st9_final[,6],st9_final[,2],st9_final[,3],st9_final[,4],st9_final[,5])
  
  
  
  
  final<-rbind(final_tiger,final_tiger2,final_tiger3,final_tiger4,final_tiger5,final_tiger6,final_tiger7)
  return(final)
}


