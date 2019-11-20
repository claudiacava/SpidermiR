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


#' @title Download mmu miRNA predicted database
#' @description SpidermiRdownload_miRNAprediction_mmu will download miRNA predicted target
#' @param mirna_list miRNA list of interest
#' @examples
#' mirna<-c('mmu-miR-708-3p')
#' list<-SpidermiRdownload_miRNAprediction_mmu(mirna_list=mirna)
#' @export
#' @import stats
#' @importFrom MAGeCKFlute TransGeneID
#' @import miRNAtap.db
#' @importFrom miRNAtap getPredictedTargets 
#' @return a dataframe with miRNA target validated interactions
SpidermiRdownload_miRNAprediction_mmu<-function(mirna_list){
  dop=list()
  dop2=list()
  for (k in  1:length(mirna_list)){
    print(paste("Processing...",mirna_list[k]))
    targets <- getPredictedTargets(mirna_list[k],species='mmu', method ='geom')
    if(is.null(targets)){
      dop2[[k]]<-mirna_list[k]
    }
    if(length(targets)!=0){
      dop[[k]]<-targets
      names(dop)[[k]]<-mirna_list[k]}
  }
  target<-dop
  lla<-list()
  for (i in 1:length(target)){
    print(paste("Mapping to gene symbol",names(target)[i]))
    s<-TransGeneID(rownames(target[[i]]), fromType =  "Entrez" , toType = "Symbol",organism = "mmu", useBiomart = FALSE,ensemblHost = "www.ensembl.org")
    a<-rep(names(target)[i],length(s))
    results<-cbind(a,toupper(s))
    lla[[i]]<-results
  }
  mat<-do.call("rbind",lla)
  
  
  return(mat)
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
    #site_mir2disease<-.url_cache$get("miRtar")
    #mir2disease<-read.delim(site_mir2disease,header = FALSE,quote = "",stringsAsFactors=FALSE)
    # querying miRNA WALK database (validated interaction miRNA-gene)
    temp <- tempfile()
    download.file(.url_cache$get("miRwalk"),temp)
    sx<-unz(temp,"hsa-vtm-gene.rdata.Rdata")
    load(sx)
    id<-t(sapply(id, '[', 1:max(sapply(id, length)))) 
    se<-int(id)
    # merging miRtar and miRNA walk information

    mir_validated_targe<-se
    site_mirtarbase<-.url_cache$get("miRTarBase")
    test <- read.xls(site_mirtarbase, quote="",stringsAsFactors=FALSE)
    pro<-as.data.frame(cbind(test$X.miRNA.,test$X.Target.Gene.))
    dem<- as.data.frame(sapply(pro, function(x) gsub("\"", "", x)))
    mir_validated_targe3<-rbind(mir_validated_targe,dem)
    mir_validated_targe4<-mir_validated_targe3[!duplicated(mir_validated_targe3), ]
    
    return(mir_validated_targe4)
}



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
  url1<-.url_cache$get("matador")
  con <- gzcon(url(url1))
  txt <- readLines(con)
  zz<-read.delim(textConnection(txt))
  zz_m<-zz[,c(5,2)]
  colnames(zz_m)<-colnames(dgib_m)
  dgib_m$gene_name<-toupper(dgib_m$gene_name)
  dgib_m$drug_name<-toupper(dgib_m$drug_name)
  zz_m$gene_name<-toupper(zz_m$gene_name)
  zz_m$drug_name<-toupper(zz_m$drug_name)
  dgidb_matador<-rbind(zz_m,dgib_m)
  dgidb_matador_tot<-dgidb_matador[!duplicated(dgidb_matador), ]
  return(dgidb_matador_tot)
}



#' @title Download drug-miRNA interactions with fisher test 
#' @description SpidermiRdownload_pharmacomir will calculate miRNA- drug interactions 
#' @param list miRNA gene target as obtained from e.g.; SpidermiRdownload_miRNAvalidate
#' @param drug parameter drug of interest
#' @examples
#' list1<-SpidermiRdownload_miRNAvalidate(validated)
#' drug="TAMOXIFEN"
#' drug_genetarget<-SpidermiRdownload_pharmacomir(list1[1:100,],drug="TAMOXIFEN")
#' @export
#' @import stats
#' @return a dataframe with miRNA target validated interactions
SpidermiRdownload_pharmacomir<-function(list,drug){
mirna_uni<-unique(list[,1])
table_pathway_enriched <- matrix(0, length(mirna_uni),5)
colnames(table_pathway_enriched) <- c("drug_target","miRNA_target","common_gene","Pvalue","FDR")
rownames(table_pathway_enriched)<- mirna_uni                     
table_pathway_enriched <- as.data.frame(table_pathway_enriched)
dgidb<-SpidermiRdownload_drug_gene(drug_gene)
a<-unique(dgidb$drug_name)
dgidb_SEL<-dgidb[dgidb$drug_name==drug,]
one_gene<-unique(dgidb_SEL$gene_name)
one_gene<-toupper(one_gene[one_gene!=""])
#i=842
#i=1
for (i in 1:length(mirna_uni)){
  print(paste0("processing...... ", mirna_uni[i], "....n. ",   i, " of......", length(mirna_uni)  ))
  mirna_name<-mirna_uni[i]
list_sel<-list[list$V1==mirna_name,]
list_sel$V2<-toupper(list_sel$V2)
genes_common_pathway_TFregulon <- as.matrix(intersect(toupper(one_gene),toupper(list_sel$V2)))
if (length(genes_common_pathway_TFregulon) != 0) {
  current_pathway_commongenes_num <- length(genes_common_pathway_TFregulon)
  allgene<-unique(dgidb$gene_name)
  seta <-  allgene %in% list_sel$V2 # tutti i deg geni
  setb <-  allgene %in% one_gene 
  ft <- fisher.test(seta,setb)
  FisherpvalueTF <- ft$p.value
  table_pathway_enriched[i,"Pvalue"] <- as.numeric(FisherpvalueTF)
  if (FisherpvalueTF < 0.05) {
     table_pathway_enriched[i,"miRNA_target"] <- nrow(list_sel)
  table_pathway_enriched[i,"drug_target"] <- length(one_gene)
  table_pathway_enriched[i,"common_gene"] <- length(genes_common_pathway_TFregulon)
  
  } 
} 
}
  table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[,"Pvalue"],decreasing = FALSE),]
  table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"Pvalue"] < 0.05 ,]
  table_pathway_enriched[,"FDR"] <- p.adjust(table_pathway_enriched[,"Pvalue"],method = "fdr")
  FDRThresh=0.05
  table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"FDR"] < FDRThresh ,]
  table_pathway_enriched <- table_pathway_enriched[order(table_pathway_enriched[,"FDR"],decreasing = FALSE),]
  n_CommonGenes=2
  table_pathway_enriched <- table_pathway_enriched[table_pathway_enriched[,"common_gene"]>n_CommonGenes,]
 return(table_pathway_enriched) 
}

 

