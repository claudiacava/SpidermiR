




#' @title Searching by biomarkers of interest with direct interaction
#' @description SpidermiRanalyze_direct_net  finds other biomarkers that are related to a set of biomarkers of interest (the input of user) with direct interations.
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param BI  a set of biomarkers of interest 
#' @export
#' @return dataframe with direct interaction of biomarkers of interest
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','FOXM1'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' biomark_of_interest<-c("hsa-let-7a","CDK","FOXO1","hsa-miR-27a")
#' GIdirect_net<-SpidermiRanalyze_direct_net(data=miRNA_cN,BI=biomark_of_interest)
SpidermiRanalyze_direct_net<-function(data,BI){
  colnames(data) <- c("gene_symbolA", "gene_symbolB")
  if(is.data.frame(data)!='TRUE'){
    data<-do.call("rbind", data)
    data<-as.data.frame(data[!duplicated(data), ]) }
  vd=list()
  vdb=list()
  j=1
  for (j in 1:length(BI)){
    x<-BI[j]
    de<-data[which(data$gene_symbolA==x),]
    deb<-data[which(data$gene_symbolB==x),]
    if(nrow(de)==0 && nrow(deb)==0){
      print(paste(BI[j],"is not in the network or please check the correct name"))}
    if(nrow(de)!=0 | nrow(deb)!=0){
      vd[[j]]<-de
      vdb[[j]]<-deb
      #print(paste(" Download genes n. ", j ," ", BI[j], " of ",length(BI), sep = ""))
    }
    }
  ds<-do.call("rbind", vd) 
  dsb<-do.call("rbind",vdb)
  dat_ex<-cbind(ds$gene_symbolA,ds$gene_symbolB)
  dat2_ex<-cbind(dsb$gene_symbolA,dsb$gene_symbolB)
  gene_inter<-as.data.frame(rbind(dat_ex,dat2_ex))
  gene_inter<-as.data.frame(gene_inter[!duplicated(gene_inter), ])
  i <- sapply(gene_inter, is.factor)
  gene_inter[i] <- lapply(gene_inter[i], as.character)
  return(gene_inter)
}






#' @title Searching by biomarkers of interest with direct interaction by ONLY the nodes in BI
#' @description SpidermiRanalyze_direct_subnetwork creates a sub network composed by ONLY the nodes in genes of interest and the edges between them
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param BI  a set of biomarkers of interest 
#' @importFrom igraph graph.data.frame induced.subgraph
#' @export
#' @return dataframe with direct interaction of biomarkers of interest
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','FOXM1'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' biomark_of_interest<-c("hsa-let-7a","CDK","FOXO1","hsa-miR-27a")
#' subnet<-SpidermiRanalyze_direct_subnetwork(data=miRNA_cN,BI=biomark_of_interest)
SpidermiRanalyze_direct_subnetwork<-function(data,BI){
  colnames(data) <- c("gene_symbolA", "gene_symbolB")
  if(is.data.frame(data)!='TRUE'){
    data<-do.call("rbind", data)
    data<-cbind(data$gene_symbolA,data$gene_symbolB)
    data<-as.data.frame(data[!duplicated(data), ]) }
  i <- sapply(data, is.factor)
  data[i] <- lapply(data[i], as.character)
  ver<-unlist(data)
  n<-unique(ver)
  s<-intersect(n,BI)
  if(length(s) !=length(BI) ){
    cv<-setdiff(BI,s)
    #BI<-setdiff(BI,cv)
    #print(paste('Error: gene,',cv, ', is not in the gene network or please check the correct name!'))
  }
  g <- graph.data.frame(data,directed=FALSE)
  g2 <- induced.subgraph(graph=g,vids=s)
  aaa<-get.data.frame(g2)
  colnames(aaa)[1] <- 'V1'
  colnames(aaa)[2] <- 'V2'
  #visualize_gi(aaa,subv)
  return(aaa)
}


#' @title Searching by biomarkers of interest and all the edges among this bunch of nodes
#' @description SpidermiRanalyze_subnetwork_neigh create a sub network composed by the nodes in BI and, if some of them are connected to other nodes (even if not in BI), take also them (include all the edges among this bunch of nodes).
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param BI  a set of biomarkers of interest 
#' @importFrom igraph graph.data.frame induced.subgraph decompose.graph get.data.frame V
#' @export
#' @return dataframe with interactions
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' biomark_of_interest<-c("hsa-let-7a","CDK","FOXO1","hsa-miR-27a")
#' GIdirect_net_neigh<-SpidermiRanalyze_subnetwork_neigh(data=miRNA_cN,BI=biomark_of_interest)
SpidermiRanalyze_subnetwork_neigh<-function(data,BI){
  data<-as.data.frame(data[!duplicated(data), ]) 
  i <- sapply(data, is.factor)
  data[i] <- lapply(data[i], as.character)
  ver<-unlist(data)
  n<-unique(ver)
  s<-intersect(n,BI)
  if(length(s) !=length(BI) ){
    cv<-setdiff(BI,s)
    print(paste('Error: gene,',cv, ', is not in the gene network or please check the correct name!'))}
  g <- graph.data.frame(data,directed=FALSE)
  sg1 <- decompose.graph(g,mode="weak")
  neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% BI)) V(s)$name else NULL})))
  g3 <- induced.subgraph(graph=g,vids=neighverts)
  aaa2<-get.data.frame(g3)
  colnames(aaa2)[1] <- 'V1'
  colnames(aaa2)[2] <- 'V2'
  #visualize2(aaa2)
  return(aaa2)}


#' @title Ranking degree centrality genes
#' @description SpidermiRanalyze_degree_centrality provides degree centrality, defined as the total number of direct neighbors for each gene.
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param cut  parameter cut is able to cut off other genes 
#' @export
#' @return dataframe with the ranked number of direct neighobors for each gene of the network
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' biomark_of_interest<-c("hsa-let-7a","CDK","FOXO1","hsa-miR-27a")
#' top10_cent<-SpidermiRanalyze_degree_centrality(miRNA_cN)
SpidermiRanalyze_degree_centrality<-function(data,cut=NULL){
  colnames(data) <- c("gene_symbolA", "gene_symbolB")
  dfer<-unlist(data)
  df<-table(dfer)
  s<-as.data.frame(df)
  if(!is.null(cut)){
    hg<-as.data.frame(s[ order(s$Freq,decreasing = TRUE ),]) 
    mol<-as.data.frame(hg[1:cut,])
    return(mol)
  }
  if(is.null(cut)){
    hg<-as.data.frame(s[ order(s$Freq,decreasing = TRUE ),]) 
    return(hg)}
}


