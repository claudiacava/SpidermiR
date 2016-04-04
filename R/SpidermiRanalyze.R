#' @title Integration of microRNA target networks.
#' @description SpidermiRanalyze_mirna_network creates a data frame with miRNA gene interaction. The user can filter the search by disease.
#' @param data  SpidermiRprepare_NET output
#' @param disease miRNA gene interaction can be filtered by disease using the parameters obatined from SpidermiRquery_disease
#' @export
#' @import stats
#' @return dataframe with miRNA gene interaction data
#' @examples
#' GS_net <- data.frame(gA=c('SMAD','MYC'),gB=c('FOXM1','KRAS'),stringsAsFactors=FALSE)
#' miRNA_NET<-SpidermiRanalyze_mirna_network(data=GS_net,disease="prostate cancer")
SpidermiRanalyze_mirna_network<-function(data,disease=NULL){
  # querying miRtar database (validated interaction miRNA-gene)
  site_mir2disease<-.url_cache$get("miRtar")
  mir2disease<-read.delim(site_mir2disease,header = FALSE,quote = "",stringsAsFactors=FALSE)
  # querying miRNA WALK database (validated interaction miRNA-gene)
  temp <- tempfile()
  download.file(.url_cache$get("miRwalk"),temp)
  sx<-unz(temp,"hsa-vtm-gene.rdata.Rdata")
  load(sx)
  id<-t(sapply(id, '[', 1:max(sapply(id, length)))) 
  se=list()
  for (j in 1:nrow(id)){
    az<-id[j,]
    de<-as.data.frame(az)
    de[,2]<-rep(rownames(id)[j],length(de))
    de<-de[complete.cases(de),]
    se[[j]]<-de
  }
  ds<-do.call("rbind", se)
  ds<-ds[c(2,1)]
  # merging miRtar and miRNA walk information
  colnames(ds) <- c("V1", "V2")
  mir2disease$V3<-NULL
  mir2disease$V4<-NULL
  mir2disease<-rbind(mir2disease,ds)
  #for the disease
  all_entries<-.url_cache$get("miR2Disease")
  disease_ref<-read.delim(all_entries,header = FALSE,quote = "",stringsAsFactors=FALSE)
  #FIND MIRNA LINK TO A PARTICULAR DISEASE
  if( !is.null(disease) ){
      AZXC<-disease_ref[disease_ref$V2==disease,]
  }
  if (is.null(disease)) {
    AZXC<-disease_ref
  }
  #rownames(AZXC)<-AZXC$V1
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  AZXC$V1<-trim(AZXC$V1)
  am<-match(mir2disease$V1,AZXC$V1)
  st3 <- cbind( mir2disease, AZXC[am,] )
  ASWQ<-na.omit(st3)
  ASWQ<-as.data.frame(ASWQ[!duplicated(ASWQ), ]) 
  vd<-list()
  if(is.data.frame(data)=='FALSE'){
  data<-do.call("rbind", data)
  }
  data<-as.data.frame(data[!duplicated(data), ]) 
  if(ncol(data)==2){
    colnames(data) <- c("gene_symbolA", "gene_symbolB")
  }
  vdb<-list()
  for (j in 1:length(ASWQ$V2)){
    x<-ASWQ$V2[j]
    de<-data[which(data$gene_symbolA==x),]
    if(nrow(de)!=0){
      de$miRNASyA<-ASWQ$V1[j]}
    vd[[j]]<-de
    deb<-data[which(data$gene_symbolB==x),]
    if(nrow(deb)!=0){
      deb$miRNASyB<-ASWQ$V1[j]}
    vdb[[j]]<-deb
    #print(paste(" Download gene n. ", j," ",ASWQ$V2[j], " and miRNA ",ASWQ$V1[j], " of ",length(ASWQ$V2), sep = ""))
  }
  ds<-do.call("rbind", vd)
  dsb<-do.call("rbind", vdb)
  dat<-cbind(ds$miRNASyA, ds$gene_symbolA)
  datb<-cbind(dsb$miRNASyB, dsb$gene_symbolB)
  dati_mirna<-as.data.frame(rbind(dat,datb))
  dati_mirna<-dati_mirna[!duplicated(dati_mirna), ]
  i <- sapply(dati_mirna, is.factor)
  dati_mirna[i] <- lapply(dati_mirna[i], as.character)
  return(dati_mirna)
}


#' @title Integration of microRNA target gene networks.
#' @description SpidermiRanalyze_mirna_gene_complnet creates a data frame with miRNA target gene interaction. The user can filter the search by disease.
#' @param data SpidermiRprepare_NET output
#' @param disease miRNA target gene interaction can be filtered by disease using the parameters obatined from SpidermiRquery_disease
#' @export
#' @import stats
#' @return dataframe with miRNA target gene interaction data
#' @examples
#' GS_net <- data.frame(gA=c('SMAD','MYC'),gB=c('FOXM1','KRAS'),stringsAsFactors=FALSE)
#' miRNA_cNT<-SpidermiRanalyze_mirna_gene_complnet(data=GS_net,disease="prostate cancer")
SpidermiRanalyze_mirna_gene_complnet<-function(data,disease=NULL){
  # querying miRtar database (validated interaction miRNA-gene)
  site_mir2disease<-.url_cache$get("miRtar")
  mir2disease<-read.delim(site_mir2disease,header = FALSE,quote = "",stringsAsFactors=FALSE)
  # querying miRNA WALK database (validated interaction miRNA-gene)
  temp <- tempfile()
  download.file(.url_cache$get("miRwalk"),temp)
  sx<-unz(temp,"hsa-vtm-gene.rdata.Rdata")
  load(sx)
  id<-t(sapply(id, '[', 1:max(sapply(id, length)))) 
  se=list()
  for (j in 1:nrow(id)){
    az<-id[j,]
    de<-as.data.frame(az)
    de[,2]<-rep(rownames(id)[j],length(de))
    de<-de[complete.cases(de),]
    se[[j]]<-de
  }
  ds<-do.call("rbind", se)
  ds<-ds[c(2,1)]
  # merging miRtar and miRNA walk information
  colnames(ds) <- c("V1", "V2")
  mir2disease$V3<-NULL
  mir2disease$V4<-NULL
  mir2disease<-rbind(mir2disease,ds)
  #for the disease
  all_entries<-.url_cache$get("miR2Disease")
  disease_ref<-read.delim(all_entries,header = FALSE,quote = "",stringsAsFactors=FALSE)
  #FIND MIRNA LINK TO A PARTICULAR DISEASE
  if(!is.null(disease)){
      AZXC<-disease_ref[disease_ref$V2==disease,]
    }
  if (is.null(disease)) {
    AZXC<-disease_ref
  }  
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  AZXC$V1<-trim(AZXC$V1)
  am<-match(mir2disease$V1,AZXC$V1)
  st3 <- cbind( mir2disease, AZXC[am,] )
  ASWQ<-na.omit(st3)
  ASWQ<-as.data.frame(ASWQ[!duplicated(ASWQ), ]) 
  list_network_mirna<-list()
  vd<-list()
  vdb<-list()
  if(is.data.frame(data)=='FALSE'){
  data<-do.call("rbind", data)
  }
  data<-as.data.frame(data[!duplicated(data), ]) 
  if(ncol(data)==2){
    colnames(data) <- c("gene_symbolA", "gene_symbolB")
  }
  for (j in 1:length(ASWQ$V2)){
    x<-ASWQ$V2[j]
    de<-data[which(data$gene_symbolA==x),]
    if(nrow(de)!=0){
      de$miRNASyA<-ASWQ$V1[j]}
    vd[[j]]<-de
    deb<-data[which(data$gene_symbolB==x),]
    if(nrow(deb)!=0){
      deb$miRNASyB<-ASWQ$V1[j]}
    vdb[[j]]<-deb
    #print(paste(" Download genes n. ", j," ",ASWQ$V2[j], " and miRNAs ",ASWQ$V1[j], " of ",length(ASWQ$V2), sep = ""))
  }
  
  ds<-do.call("rbind", vd)
  dsb<-do.call("rbind", vdb)
  dat<-cbind(ds$miRNASyA,ds$gene_symbolA)
  dat2<-cbind(ds$gene_symbolA,ds$gene_symbolB)
  zzz<-rbind(dat,dat2)
  datb<-cbind(dsb$miRNASyB,dsb$gene_symbolB)
  dat2b<-cbind(dsb$gene_symbolA,dsb$gene_symbolB)
  cv<-rbind(datb,dat2b)
  fb<-as.data.frame(rbind(zzz,cv))
  s_n<-as.data.frame(fb[!duplicated(fb), ])
  i <- sapply(s_n, is.factor)
  s_n[i] <- lapply(s_n[i], as.character)
  return(s_n)
}




#' @title Integration of Extracellular/Circulating miRNA
#' @description SpidermiRanalyze_mirna_extra_cir creates a data frame with miRNA target or miRNA target gene interaction.
#' @param data SpidermiRanalyze_mirna_gene_complnet output or network with miRNA in the first column of the dataframe
#' @param type The user using the following parameteres can specify the network type
#' \tabular{ll}{
#' mT \tab   to obtain a microRNA target interactions \cr
#' mCT  \tab   to obtain a microRNA target gene complete iteractions\cr}
#' @export
#' @import stats
#' @return dataframe with Extracellular/Circulating miRNA target interaction data
#' @examples
#'miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#'miRNA_NET_ext_circmT<-SpidermiRanalyze_mirna_extra_cir(data=miRNA_cN,"mT")
SpidermiRanalyze_mirna_extra_cir<-function(data,type=NULL){
  site<-.url_cache$get("mirandola")
  mirandola<-read.delim(site,header = TRUE,quote = "",stringsAsFactors=FALSE)
  colnames(data) <- c("V1", "V2")
 vd=list()
  for (j in 1:length(mirandola$RNA_name)){
    x<-mirandola$RNA_name[j]
    de<-data[which(data$V1==x),]
    if(nrow(de)!=0){
      de$miRNASyA<-mirandola$RNA_name[j]
    vd[[j]]<-de}}
    ds<-do.call("rbind", vd)
    data2<-as.data.frame(ds[!duplicated(ds), ]) 
  ASWQ<-as.data.frame(lapply(data2, as.character),stringsAsFactors=FALSE)
  cd<-as.data.frame(cbind(ASWQ$V1,ASWQ$V2))
  cd<-cd[!duplicated(cd), ]
  if(!is.null(type)){
    if (type=="mT"){
      i <- sapply(cd, is.factor)
      cd[i] <- lapply(cd[i], as.character)
      return(cd)}
    if (type=="mCT"){
      nh<-data[-grep("hsa",data$V1),]
      def<-rbind(nh,cd)
      i <- sapply(def, is.factor)
      def[i] <- lapply(def[i], as.character)
      return(def)}
  }
  if (is.null(type)) {
    nh<-data[-grep("hsa",data$V1),]
    def<-rbind(nh,cd)
    i <- sapply(def, is.factor)
    def[i] <- lapply(def[i], as.character)
    return(def)}
}


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


#' @title Find community detection
#' @description SpidermiRanalyze_Community_detection try to find dense subgraphs in directed or undirected graphs, by optimizing some criteria.
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param type  with the parameter type the user can choose the algorithm  to calculate the community structure
#' EB edge.betweenness.community
#' FC fastgreedy.community
#' WC walktrap.community
#' SC spinglass.community
#' LE leading.eigenvector.community
#' LP label.propagation.community
#' @importFrom igraph graph.data.frame 
#' @importFrom igraph induced.subgraph 
#' @importFrom igraph  decompose.graph 
#' @importFrom igraph get.data.frame 
#' @importFrom igraph delete.edges
#' @importFrom igraph edge.betweenness.community 
#' @importFrom igraph fastgreedy.community 
#' @importFrom igraph walktrap.community 
#' @importFrom igraph spinglass.community 
#' @importFrom igraph leading.eigenvector.community 
#' @importFrom igraph label.propagation.community
#' @importFrom igraph ecount
#' @importFrom igraph clusters
#' @importFrom igraph modularity
#' @export
#' @return a list of clusters with their number of genes
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' comm<-  SpidermiRanalyze_Community_detection(data=miRNA_cN,type="FC") 
SpidermiRanalyze_Community_detection<-function(data,type){
  colnames(data) <- c("gene_symbolA", "gene_symbolB")
  g <- graph.data.frame(data,directed=FALSE)
  if (type=="EB"){
    eb <- edge.betweenness.community(g,directed=FALSE)
  }
  if (type=="FC"){
    eb <- fastgreedy.community(g)}
  if (type=="WC"){
    eb <- walktrap.community(g)
  }
  if (type=="SC"){
    eb <- spinglass.community(g)
  }
  if (type=="LE"){
    eb <- leading.eigenvector.community(g)
  }
  if (type=="LP"){
    eb <- label.propagation.community(g)
  }
  mods <- sapply(0:ecount(g), function(i){
    g2 <- delete.edges(g, eb$removed.edges[seq(length=i)])
    cl <- clusters(g2)$membership
    modularity(g,cl)})
  g2<-delete.edges(g, eb$removed.edges[seq(length=which.max(mods)-1)])
  #To determine the number of vertices in each cluster
  print(clusters(g2)$csize)
  return(g2)}





#' @title Find the network of community detection
#' @description SpidermiRanalyze_Community_detection_net find the direct interactions from a specific community
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param comm_det SpidermiRanalyze_Community_detection
#' @param size the index of community detection obtained from SpidermiRanalyze_Community_detection
#' @export
#' @return dataframe with the interatcions
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' comm<-  SpidermiRanalyze_Community_detection(data=miRNA_cN,type="FC") 
#' cd_net<-SpidermiRanalyze_Community_detection_subnet(data=miRNA_cN,comm_det=comm,size=1)
SpidermiRanalyze_Community_detection_subnet<-function(data,comm_det,size){
  colnames(data) <- c("gene_symbolA", "gene_symbolB")
  #want to find the vertices in cluster 
  if(is.data.frame(comm_det)!='TRUE'){
    z<-as.data.frame(which(clusters(comm_det)$membership == size))
    a<-rownames(z)
    as<-SpidermiRanalyze_direct_subnetwork(data,a)
    return(as)
  }}



#' @title Find the network of community detection and direct biormarker
#' @description SpidermiRanalyze_direct_net find the direct interactions from a specific community
#' @param data  SpidermiRanalyze_mirna_network output or SpidermiRanalyze_mirna_gene_complnet
#' @param comm_det SpidermiRanalyze_Community_detection
#' @param size the index of community detection obtained from SpidermiRanalyze_Community_detection
#' @export
#' @return dataframe with the interatcions
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' comm<-  SpidermiRanalyze_Community_detection(data=miRNA_cN,type="FC") 
#' cd_net<-SpidermiRanalyze_Community_detection_net(data=miRNA_cN,comm_det=comm,size=1)
SpidermiRanalyze_Community_detection_net<-function(data,comm_det,size){
  colnames(data) <- c("gene_symbolA", "gene_symbolB")
  #want to find the vertices in cluster 
  if(is.data.frame(comm_det)!='TRUE'){
    z<-as.data.frame(which(clusters(comm_det)$membership == size))
    a<-rownames(z)
    as<-SpidermiRanalyze_direct_net(data,a)
    return(as)
  }}






#' @title Community detection from biomarkers of interest
#' @description SpidermiRanalyze_Community_detection_bi find the cluster with biomarkers of interest
#' @param data SpidermiRanalyze_Community_detection output
#' @param BI a set of biomarkers of interest
#' @export
#' @return a list with the cluster for each biomarkers of interest
#' @examples
#' miRNA_cN <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' comm<-  SpidermiRanalyze_Community_detection(data=miRNA_cN,type="FC") 
#' biomark_of_interest<-c("hsa-let-7a","CDK","FOXO1","hsa-miR-27a")
#' mol<-SpidermiRanalyze_Community_detection_bi(data=comm,BI=biomark_of_interest)
SpidermiRanalyze_Community_detection_bi<-function(data,BI){
  bg=list()
  for (j in 1:length(clusters(data)$csize)){ 
    z<-as.data.frame(which(clusters(data)$membership == j))
    z[1]<-rownames(z)
    z[2]<-j
    bg[[j]]<-z}
  ds<-do.call("rbind", bg) 
  ge=list()
  for (j in 1:length(BI)){
    x<-BI[j]
    de<-ds[which(rownames(ds)==x),]
    ge[[j]]<-de    
    if(nrow(de)!=0){
      print(paste(x,"in the cluster n",de$V2))}
    if(nrow(de)==0){
      print(paste(x,"doesn't find in the cluster"))}
  }
  dsn<-do.call("rbind", ge)
  return(dsn)
} 


#' @title Integration of pharmacomiR in the network
#' @description SpidermiRanalyze_mirnanet_pharm integrates both miRNA targeting of the gene and the gene-drug interaction from PharmacomiR database in the network
#' @param mir_ph SpidermiRdownload_pharmacomir output
#' @param net a network data (e.g. SpidermiRanalyze_mirna_network or SpidermiRanalyze_mirna_gene_complnet output)
#' @export
#' @return a dataframe with the integation of network and pharmacomiR data
#' @examples
#' mir_p <-data.frame(gA=c('hsa-let-7a','CASP3'),gB=c('CASP3','paclitaxel'),stringsAsFactors=FALSE)
#' net_p <-data.frame(gA=c('hsa-let-7a','hsa-miR-300'),gB=c('FOXM1','KPNA4'),stringsAsFactors=FALSE)
#' mol<-SpidermiRanalyze_mirnanet_pharm(mir_ph=mir_p,net=net_p)
SpidermiRanalyze_mirnanet_pharm<-function(mir_ph,net){
  sss<-rbind(net,mir_ph)
  s2<-as.data.frame(sss[!duplicated(sss), ]) 
  return(s2)
}
