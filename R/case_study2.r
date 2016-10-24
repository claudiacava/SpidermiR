#' @title Loading 1 network of Protein Interactions (PI)
#' @description Case_Study2_loading_1_network loads PI in HomoSapiens
#' @param species  variable
#' @export
#' @return dataframe with interactions
#' @examples
#' \dontrun{
#' a2<-Case_Study2_loading_1_network(species)}
Case_Study2_loading_1_network<-function(species){
org<-SpidermiRquery_species(species)
net_PHint<-SpidermiRquery_spec_networks(organismID = org[6,],network = "PHint")
out_net<-SpidermiRdownload_net(net_PHint)
geneSymb_net<-SpidermiRprepare_NET(organismID = org[6,],data = out_net)
ds<-do.call("rbind", geneSymb_net)
data1<-as.data.frame(ds[!duplicated(ds), ]) 
sdas<-cbind(data1$gene_symbolA,data1$gene_symbolB)
sdas<-as.data.frame(sdas[!duplicated(sdas), ]) 
m<-c(data1$gene_symbolA)
m2<-c(data1$gene_symbolB)
s<-c(m,m2)
fr<- unique(s)
network="PHint"
print(paste("Downloading of 1 ",network, " network ","in ",org[6,]," with number of nodes: ",length(fr)," and number of edges: ",nrow(sdas), sep = ""))
return(geneSymb_net)
}



#' @title Loading 2 network of Protein Interactions (PI) with miRNAs
#' @description Case_Study2_loading_2_network loads PI in HomoSapiens with miRNAs already found as deregulated in BC (only interaction miRNA-gene)
#' @param data  output of Case_Study2_loading_1_network
#' @export
#' @return dataframe with interactions
#' @examples
#' \dontrun{
#' b2<-Case_Study2_loading_2_network(data=a2)
#' }
Case_Study2_loading_2_network<-function(data){
miRNA_NET<-SpidermiRanalyze_mirna_network(data,disease="breast cancer",miR_trg="val")
m2<-c(miRNA_NET$V1)
m3<-c(miRNA_NET$V2)
s2<-c(m2,m3)
fr2<- as.data.frame(unique(s2))
print(paste("Downloading of 2 network with the integration of miRNA-gene interaction with number of nodes ", nrow(fr2)," and number of edges ", nrow(miRNA_NET), sep = ""))
return(miRNA_NET)
}
#a<-fr2[grep("hsa",fr2$`unique(s2)`),]




#' @title Loading 3 network of Protein Interactions (PI) with miRNAs
#' @description Case_Study2_loading_2_network loads PI in HomoSapiens with miRNAs already found as deregulated in BC (only interaction miRNA-gene)
#' @param sdas  output of Case_Study2_loading_1_network
#' @param miRNA_NET  output of Case_Study2_loading_2_network
#' @export
#' @return dataframe with interactions
#' @examples
#' \dontrun{
#' c2<-Case_Study2_loading_3_network(sdas=a2,miRNA_NET=b2)
#' }
Case_Study2_loading_3_network<-function(sdas,miRNA_NET){
ds<-do.call("rbind", sdas)
  data1<-as.data.frame(ds[!duplicated(ds), ]) 
  sdas<-cbind(data1$gene_symbolA,data1$gene_symbolB)
  sdas<-as.data.frame(sdas[!duplicated(sdas), ]) 
topwhol<-SpidermiRanalyze_degree_centrality(sdas)
topwhol_mirna<-SpidermiRanalyze_degree_centrality(miRNA_NET)

miRNA_degree<-topwhol_mirna[grep("hsa",topwhol_mirna$dfer),]
seq_gd<-as.data.frame(seq(1, 15400, by = 50))
even<-seq_gd[c(F,T),]
even2<-even
odd<-seq_gd[c(T,F),]
odd2<-odd[-1]
odd2[154]<-15400
f<-cbind(even2,odd2-1)

SQ<-cbind(odd,even-1)

h<-as.data.frame(rbind(f,SQ))
SQ <- as.data.frame(h[order(h$even2,decreasing=FALSE),])

table_pathway_enriched <- matrix(1, nrow(SQ),4)
colnames(table_pathway_enriched) <- c("interval min","interval max","gene","miRNA")
table_pathway_enriched <- as.data.frame(table_pathway_enriched)

j=1
for (j in 1:nrow(SQ)){ 
  a<-SQ$even2[j]
  b<-SQ$V2[j]
  d<-c(a,b)
gene_degree10<-topwhol[a:b,]
vfg<-rbind(miRNA_degree[1:10,],gene_degree10)
subnet<-SpidermiRanalyze_direct_subnetwork(data=miRNA_NET,BI=vfg$dfer)

table_pathway_enriched[j,"interval min"] <- d[1]
table_pathway_enriched[j,"interval max"] <- d[2]
s<-unique(subnet$V1)
x<-unique(subnet$V2)
table_pathway_enriched[j,"miRNA"]<-length(s)
table_pathway_enriched[j,"gene"]<-length(x)
}

df<-cbind(table_pathway_enriched$gene,table_pathway_enriched$miRNA)
rownames(df)<-table_pathway_enriched$`interval max`
categories <- c("protein", "miRNA")
colors <- c("green", "magenta")
op <- par(mar = c(5, 5, 4, 2) + 0.1)
matplot(df, type="l",col=colors,xlab = "N of Clusters",main = "",ylab = "Interactions",cex.axis=2,cex.lab=2,cex.main=2)
legend("topright", col=colors, categories, bg="white", lwd=1,cex=2)
j=1
a<-SQ$even2[j]
b<-SQ$V2[j]
d<-c(a,b)
gene_degree10<-topwhol[a:b,]
vfg<-rbind(miRNA_degree[1:10,],gene_degree10)
subnet<-SpidermiRanalyze_direct_subnetwork(data=miRNA_NET,BI=vfg$dfer)
m2<-c(subnet$V1)
m3<-c(subnet$V2)
s2<-c(m2,m3)
fr2<- as.data.frame(unique(s2))


print(paste("Downloading of 3 network with proteins and miRNAs with highest degree centrality with  ", nrow(fr2)," nodes and number of edges ", nrow(subnet), sep = ""))


return(subnet)
}





#SpidermiRvisualize_3Dbarplot<-function(Edges_1net,Edges_2net,Edges_3net,NODES_1net,NODES_2net,NODES_3net,nmiRNAs_1net,nmiRNAs_2net,nmiRNAs_3net,ng_1net,ng_2net,ng_3net){
 # y <- as.factor(c("1net","2net","3net","1net","2net","3net","1net","2net","3net","1net","2net","3net"))
  #x<-c(round(log(Edges_1net)),round(log(Edges_2net)),round(log(Edges_3net)),round(log(NODES_1net)),round(log(NODES_2net)),round(log(NODES_3net)),round(log(nmiRNAs_1net)),round(log(nmiRNAs_2net)),round(log(nmiRNAs_3net)),round(log(ng_1net)),round(log(ng_2net)),round(log(ng_3net)))
  #z<-as.factor(c("edges","edges","edges","nodes","nodes","nodes","miRNAs","miRNAs","miRNAs","proteins","proteins","proteins"))
  #cloud(x
   #     ~y+z, panel.3d.cloud=panel.3dbars, col.facet='green',  screen = list(z = -50, x = -30),
    #    xbase=0.8, ybase=0.3, scales=list(arrows=FALSE, col=1), 
     #   par.settings = list(axis.line = list(col = "transparent")))
#}



#SpidermiRvisualize_3Dbarplot(Edges_1net=18903,Edges_2net=1001,Edges_3net=12,NODES_1net=15407,NODES_2net=830,NODES_3net=15,nmiRNAs_1net=0,nmiRNAs_2net=62,nmiRNAs_3net=7,ng_1net=15407,ng_2net=768,ng_3net=8)




