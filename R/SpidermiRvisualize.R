#' @title Visualize results obtained by SpidermiR analysis
#' @description Visualize the network
#' @param data The input data is a dataframe containing network data.
#' @examples 
#' miRNA_cNET <-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' SpidermiRvisualize_mirnanet(data=miRNA_cNET)
#' @importFrom networkD3 forceNetwork
#' @importFrom networkD3 simpleNetwork
#' @importFrom networkD3 JS
#' @export
#' @return 3D graphic
SpidermiRvisualize_mirnanet<-function(data){
  colnames(data) <- c("V1", "V2")
  if( length(grep("hsa",data$V1)) !=0 | length(grep("hsa",data$V2)) !=0){
    IDs2 <- sort.int(unique(c(data$V1, data$V2)),decreasing=FALSE)
    IDs2 <- data.frame(ID= seq_along(IDs2)-1, name= IDs2)
    dataIDs2 <- merge(data, IDs2, by.x= "V1", by.y= "name")
    dataIDs2$V1 <- dataIDs2$ID  
    dataIDs2$ID <- NULL 
    dataIDs2 <- merge(dataIDs2, IDs2, by.x= "V2", by.y= "name")
    dataIDs2$V2 <- dataIDs2$ID
    dataIDs2$ID <- NULL
    dataIDs2 <- dataIDs2[,c("V1", "V2")]
    data2 <- data[order(data$V1, data$V2,decreasing=FALSE),]
    dataIDs2 <- dataIDs2[order(dataIDs2$V1, dataIDs2$V2,decreasing=FALSE),]
    att<-as.data.frame(sort(unique(unlist(data)),decreasing=FALSE))
    colnames(att)[1]<-"v1"
    att$v2 <- "" 
    att$v2<-replace(att$v2, att$v2 == "", "gene")
    att$v2[(grep("[a-z]",att$v1))]<- "Pharmaco"
    att$v2[grep("hsa" ,att$v1)]<- "miRNA"
    att$v2[grep("orf" ,att$v1)]<- "gene"
    
    #att$v3<-NULL
    att<-as.data.frame(att[order(att$v2,decreasing = FALSE ),]) 
    i <- sapply(att, is.factor)
    att[i] <- lapply(att[i], as.character)
    colnames(att)[1]<-"name"
    colnames(att)[2]<-"Group"
    attr2 <- merge(att, IDs2)
    # Order rows
    attr2 <- attr2[order(attr2$ID),]
    #attr2<-as.data.frame(attr2[ order(attr2$Group,decreasing = FALSE ),]) 
    return(
      forceNetwork(Links = dataIDs2, Nodes = attr2, Source = "V1", Target = "V2", NodeID = "name", Group= "Group",height = 
                     800, width = 800, opacity = 1, zoom = FALSE, bounded = TRUE, legend= TRUE, opacityNoHover= 0.5,
                   colourScale=JS("d3.scale.category10()"),fontSize = 12)
    )
  }
  if( length(grep("hsa",data$V1)) ==0){
    .SpidermiRvisualize_gene(data)
  }
}


#' @title Visualize results obtained by SpidermiR analysis starting form a set of biomarker of interest
#' @description Visualize miRNA-target interaction and miRNA-target-gene starting from a set of biomarker of interest
#' @param data The input data is a dataframe containing network data.
#' @param BI a set of biomarkers of interest
#' @examples 
#' miRNA_cNET <-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' biomark_of_interest<-c("hsa-let-7a","CDK","FOXO1","hsa-miR-27a")
#' SpidermiRvisualize_BI(data=miRNA_cNET,BI=biomark_of_interest)
#' @importFrom networkD3 forceNetwork 
#' @importFrom networkD3 simpleNetwork
#' @importFrom networkD3 JS
#' @export
#' @return 3D graphic

SpidermiRvisualize_BI<-function(data,BI){
  colnames(data)[1]<-"V1"
  colnames(data)[2]<-"V2"
  IDs2 <- sort(unique(c(data$V1, data$V2)))
  IDs2 <- data.frame(ID= seq_along(IDs2)-1, name= IDs2)
  dataIDs2 <- merge(data, IDs2, by.x= "V1", by.y= "name")
  dataIDs2$V1 <- dataIDs2$ID  
  dataIDs2$ID <- NULL 
  dataIDs2 <- merge(dataIDs2, IDs2, by.x= "V2", by.y= "name")
  dataIDs2$V2 <- dataIDs2$ID
  dataIDs2$ID <- NULL
  dataIDs2 <- dataIDs2[,c("V1", "V2")]
  data2 <- data[order(data$V1, data$V2),]
  dataIDs2 <- dataIDs2[order(dataIDs2$V1, dataIDs2$V2),]
  att<-as.data.frame(sort(unique(unlist(data))))
  colnames(att)[1]<-"v1"
  att$v2 <- "" 
  for (j in 1:length(BI)){
    att$v2[grep(BI[j] ,att$v1)]<- "biomarker of interest"}
  att$v2<-replace(att$v2, att$v2 == "", "biomarker interaction")
  #att$v1 <- "name"
  i <- sapply(att, is.factor)
  att[i] <- lapply(att[i], as.character)
  colnames(att)[1]<-"name"
  colnames(att)[2]<-"Group"
  attr2 <- merge(att, IDs2)
  # Order rows
  attr2 <- attr2[order(attr2$ID),]
  return(
    forceNetwork(Links = dataIDs2, Nodes = attr2, Source = "V1", Target = "V2", NodeID = "name", Group= "Group",height = 800, width = 800, opacity = 1, zoom = FALSE, bounded = TRUE, legend= TRUE, opacityNoHover= 0.5,
                 colourScale=JS("d3.scale.category20c()")
                 ,fontSize = 12))
}



#' @title Visualize results obtained by SpidermiRanalyze_mirna_network 
#' @description It shows a plot with miRNAs and the number of their targets in the network
#' @param miRnet The input data is a dataframe containing miRNA network data (e.g. output of SpidermiRanalyze_mirna_network.
#' @examples 
#' cd<-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' SpidermiRvisualize_plot_target(miRnet=cd)
#' @importFrom ggplot2 ggplot 
#' @importFrom ggplot2 aes 
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 theme
#' @importFrom gridExtra grid.arrange
#' @export
#' @return plot
SpidermiRvisualize_plot_target<-function(miRnet){
p<-table(miRnet[,1])
as<-as.data.frame(p)
D<-as[order(as$Freq,decreasing=TRUE),]
names(D)[1]<-"miRNAs"
names(D)[2]<-"miRNAs_target"
D$miRNAs<-factor(D$miRNAs, levels=D[order(D$miRNAs_target),"miRNAs"])


return(grid.arrange(ggplot(D, aes(x = miRNAs, y = miRNAs_target)) + geom_bar(stat = "identity")+ coord_flip(),ncol=2) + theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x  = element_text(angle=360, vjust=0.5, size=16),axis.title.y = element_text(face="bold", size=16),
          axis.text.y  = element_text(angle=360, vjust=0.5, size=16)))
}
