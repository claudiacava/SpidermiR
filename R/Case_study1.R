#' @title Loading 1 network of shared protein domanin
#' @description Case_Study1_loading_1_network loads shared protein domain in HomoSapiens
#' @param species  variable
#' @export
#' @return dataframe with interactions
#' @examples
#' \dontrun{
#' a<-Case_Study1_loading_1_network(species)
#' }
Case_Study1_loading_1_network<-function(species){
org<-SpidermiRquery_species(species)
net_shar_prot<-SpidermiRquery_spec_networks(organismID = org[6,],
                                            network = "SHpd")
out_net<-SpidermiRdownload_net(net_shar_prot)
geneSymb_net<-SpidermiRprepare_NET(organismID = org[6,],data = out_net)
ds<-do.call("rbind", geneSymb_net)
data2<-as.data.frame(ds[!duplicated(ds), ]) 
m<-c(data2$gene_symbolA)
m2<-c(data2$gene_symbolB)
s<-c(m,m2)
fr<- unique(s)
network = "SHpd"
print(paste("Downloading of 1 ",network, " network ","in ",org[6,]," with number of nodes: ",length(fr)," and number of edges: ",nrow(data2), sep = ""))
return(geneSymb_net)
}







#' @title Loading 2 network of shared protein domanin
#' @description Case_Study1_loading_2_network loads shared protein domain in HomoSapiens including only miRNAs already found deregulated in PC,
#' @param data  the output of Case_Study1_loading_1_network
#' @export
#' @return dataframe with selected interactions
#' @examples
#' #' \dontrun{
#' b<-Case_Study1_loading_2_network(data=a)
#' }
Case_Study1_loading_2_network<-function(data){
miRNA_complNET<-SpidermiRanalyze_mirna_gene_complnet(data,disease="prostate cancer",miR_trg="val")
m2<-c(miRNA_complNET$V1)
m3<-c(miRNA_complNET$V2)
s2<-c(m2,m3)
fr2<- as.data.frame(unique(s2))
print(paste("Downloading of 2 network with the integration of miRNA-gene-gene interaction with number of nodes ", nrow(fr2)," and number of edges ", nrow(miRNA_complNET), sep = ""))
return(miRNA_complNET)
}





#' @title Loading 3 network of shared protein domanin
#' @description Case_Study1_loading_3_network loads shared protein domain in HomoSapiens including only the DEGs with a direct interaction among them
#' @param data  the output of Case_Study1_loading_2_network
#' @param dataFilt  TCGA matrix
#' @param dataClin  clinical data matrix
#' @importFrom TCGAbiolinks TCGAquery_SampleTypes 
#' @export
#' @return dataframe with selected interactions
#' @examples
#' \dontrun{
#' c<-Case_Study1_loading_3_network(data=b,dataFilt=dataFilt,dataClin=dataClin)
#' }
Case_Study1_loading_3_network<-function(data,dataFilt,dataClin){
#dataClin<-read.delim("C:/Users/UserInLab05/Google Drive/MIRNA AND GENEMANIA/1 case study/nationwidechildrens.org_PRAD.bio.Level_2.0.54.0/nationwidechildrens.org_clinical_patient_prad.txt")
#load("C:/Users/UserInLab05/SpidermiR_release/data/dataFilt.rda")
#take only gleason > o = 7
highstage <- dataClin[grep("7|8|9|10", dataClin$gleason_score), ]
# consider only barcode and stage
highstage<-highstage[,c("bcr_patient_barcode","gleason_score")]
highstage<-t(highstage)
samples_hight<-highstage[1,2:ncol(highstage)]
dataSmTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                  typesample = "TP")
# Which samples are Solid Tissue Normal
dataSmNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                  typesample ="NT")
colnames(dataFilt)<-substr(colnames(dataFilt),1,12)
#replace barcode of dataSmNT withthe same of clinical data
se<-substr(dataSmTP, 1, 12)
#find common patients between samples_highstage e miRNAdata
common<-intersect(colnames(dataFilt),samples_hight)
dataSmNT<-substr(dataSmNT, 1, 12)
sub_net2<-SpidermiRanalyze_DEnetworkTCGA(data,TCGAmatrix=dataFilt,tumour=common,normal=dataSmNT)
#geneSymb_net2<-cbind(data2$gene_symbolA,data2$gene_symbolB)
#sub_net<-SpidermiRanalyze_direct_subnetwork(data=b,BI=rownames(dataDEGs))
ft<-sub_net2$V1
ft1<-sub_net2$V2
fgt<-c(ft,ft1)
miRNA_NET<-SpidermiRanalyze_mirna_network(sub_net2,disease="prostate cancer",miR_trg="val")
TERZA_NET<-rbind(miRNA_NET,sub_net2)
print(paste("In the 3 network we found",length(unique(miRNA_NET$V1)), " miRNAs and ", length(unique(fgt)), " genes with ", nrow(TERZA_NET), " edges " ))
return(TERZA_NET)
}



#' @title Loading 4 network of shared protein domanin
#' @description Case_Study1_loading_4_network loads network community with the higher number of elements in the 3 network
#' @param TERZA_NET  the output of Case_Study1_loading_3_network
#' @export
#' @return dataframe with selected interactions
#' @examples
#' #' \dontrun{
#' d<-Case_Study1_loading_4_network(TERZA_NET=c)
#' }
Case_Study1_loading_4_network<-function(TERZA_NET){
comm<-  SpidermiRanalyze_Community_detection(data=TERZA_NET,type="FC")
#SpidermiRvisualize_mirnanet(TERZA_NET)
cd_net<-SpidermiRanalyze_Community_detection_net(data=TERZA_NET,comm_det=comm,size=5)
ft<-cd_net$V1
ft1<-cd_net$V2
fgt<-c(ft,ft1)
print(paste("In the 4 network we found",length(unique(fgt)), " nodes and ", nrow(cd_net), " edges " ))
return(cd_net)
}


#SpidermiRvisualize_direction(data=d)
#SpidermiRvisualize_mirnanet(cd_net)
#SpidermiRvisualize_adj_matrix(data=cd_net)
#SpidermiRvisualize_degree_dist(data=cd_net)




#SpidermiRvisualize_3Dbarplot<-function(Edges_1net,Edges_2net,Edges_3net,Edges_4net,NODES_1net,NODES_2net,NODES_3net,NODES_4net,nmiRNAs_1net,nmiRNAs_2net,nmiRNAs_3net,nmiRNAs_4net,nps_1net,nps_2net,nps_3net,nps_4net){
 # y <- as.factor(c("1net","2net","3net","4net","1net","2net","3net","4net","1net","2net","3net","4net","1net","2net","3net","4net"))
  #x<-c(round(log(Edges_1net)),round(log(Edges_2net)),round(log(Edges_3net)),round(log(Edges_4net)),round(log(NODES_1net)),round(log(NODES_2net)),round(log(NODES_3net)),round(log(NODES_4net)),round(log(nmiRNAs_1net)),round(log(nmiRNAs_2net)),round(log(nmiRNAs_3net)),round(log(nmiRNAs_4net)),round(log(nps_1net)),round(log(nps_2net)),round(log(nps_3net)),round(log(nps_4net)))
  #z<-as.factor(c("edges","edges","edges","edges","nodes","nodes","nodes","nodes","miRNAs","miRNAs","miRNAs","miRNAs","proteins","proteins","proteins","proteins"))
  #cloud(x
  #      ~y+z, panel.3d.cloud=panel.3dbars, col.facet='green',  screen = list(z = -50, x = -30),
   #     xbase=0.8, ybase=0.3, scales=list(arrows=FALSE, col=1), 
   #     par.settings = list(axis.line = list(col = "transparent")))
#}



#SpidermiRvisualize_3Dbarplot(Edges_1net=1041003,Edges_2net=100016,Edges_3net=86,Edges_4net=59,NODES_1net=16502,NODES_2net=13338,NODES_3net=80,NODES_4net=47,nmiRNAs_1net=0,nmiRNAs_2net=74,nmiRNAs_3net=13,nmiRNAs_4net=5,nps_1net=16502,nps_2net=13264,nps_3net=67,nps_4net=42)






