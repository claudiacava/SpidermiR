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
    ColourScale <- 'd3.scale.ordinal()
            .domain(["gene", "Pharmaco","miRNA"])
  .range(["#0096ff", "#00b34a","#ff6900"]);'
    return(
      forceNetwork(Links = dataIDs2, Nodes = attr2, Source = "V1", Target = "V2", NodeID = "name", Group= "Group",height = 
                     2000, width = 2000, opacity = 1, zoom = TRUE, bounded = TRUE, legend= TRUE, opacityNoHover= 0.5,
                   colourScale=JS(ColourScale),fontSize = 16)
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
  att$v2[(grep("[a-z]",att$v1))]<- "Pharmaco"
  att$v2[grep("hsa" ,att$v1)]<- "miRNA"
  
  for (j in 1:length(BI)){
    att$v2[grep(BI[j] ,att$v1)]<- "biomarker of interest"}

  att$v2[grep("orf" ,att$v1)]<- "gene"
  att$v2<-replace(att$v2, att$v2 == "", "gene")
  
  
  
  #att$v1 <- "name"
  i <- sapply(att, is.factor)
  att[i] <- lapply(att[i], as.character)
  colnames(att)[1]<-"name"
  colnames(att)[2]<-"Group"
  attr2 <- merge(att, IDs2)
  # Order rows
  attr2 <- attr2[order(attr2$ID),]
  ColourScale <- 'd3.scale.ordinal()
.domain(["gene", "Pharmaco","miRNA","biomarker of interest"])
  .range(["#0096ff", "#00b34a","#ff6900","#ffe900"]);'

  return(
    forceNetwork(Links = dataIDs2, Nodes = attr2, Source = "V1", Target = "V2", NodeID = "name", Group= "Group",height = 
                   1000, width = 1000, opacity = 1, zoom = FALSE, bounded = TRUE, legend= TRUE, opacityNoHover= 0.5,
                 colourScale=JS(ColourScale),fontSize = 16))
}



#' @title Visualize results obtained by SpidermiRanalyze_mirna_network 
#' @description It shows a plot with miRNAs and the number of their targets in the network
#' @param data The input data is a dataframe containing miRNA network data (e.g. output of SpidermiRanalyze_mirna_network.
#' @examples 
#' cd<-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' SpidermiRvisualize_plot_target(data=cd)
#' @importFrom ggplot2 ggplot 
#' @importFrom ggplot2 aes 
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 theme
#' @importFrom gridExtra grid.arrange
#' @export
#' @return plot
SpidermiRvisualize_plot_target<-function(data){
  a<-c(data[,1],data[,2])

p<-table(a)

as<-as.data.frame(p)



D<-as[order(as$Freq,decreasing=TRUE),]
names(D)[1]<-"miRNAs"
names(D)[2]<-"mRNA_target"
D$miRNAs<-factor(D$miRNAs, levels=D[order(D$mRNA_target),"miRNAs"])


return(grid.arrange(ggplot(D, aes(x = miRNAs, y = mRNA_target)) + geom_bar(stat = "identity")+ coord_flip(),ncol=2) + theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x  = element_text(angle=360, vjust=0.5, size=16),axis.title.y = element_text(face="bold", size=16),
          axis.text.y  = element_text(angle=360, vjust=0.5, size=16)))
}



#' @title plots the degree distribution of the network
#' @description It shows a plot of the degree distribution of the network
#' @param data The input data is a network
#' @examples 
#' cd<-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' SpidermiRvisualize_degree_dist(data=cd)
#' @importFrom igraph graph.data.frame degree.distribution 
#' @export
#' @return plot
SpidermiRvisualize_degree_dist<-function(data){
  g <- graph.data.frame(data,directed=FALSE)
dd <- degree.distribution(g, cumulative=T, mode="all")
return(plot(dd, pch=19, cex=1, col="orange", xlab="Degree", ylab="Cumulative Frequency")
)
}


#' @title plots the adjacency matrix of the network
#' @description It shows a plot OF the adjacency matrix of the network
#' @param data The input data is a network
#' @examples 
#' cd<-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' SpidermiRvisualize_adj_matrix(data=cd)
#' @importFrom igraph graph.data.frame as_adjacency_matrix ecount get.adjacency E E<-
#' @importFrom grDevices colorRampPalette
#' @importFrom stats  runif
#' @importFrom gplots heatmap.2
#' @export
#' @return plot
SpidermiRvisualize_adj_matrix<-function(data){
  #Convert a graph to an adjacency matrix
  g <- graph.data.frame(data,directed=FALSE)
  as_adjacency_matrix(g)
  E(g)$weight <- runif(ecount(g))
  netm <- get.adjacency(g, attr="weight", sparse=F)
  scaleyellowred <- colorRampPalette(c("red","yellow"), space = "rgb")(100)
  lhei <- c(7, 4)
  lwid=c(7, 4, 4 )
  return(heatmap.2(netm,dendrogram="none",Rowv=FALSE, Colv=FALSE,
            col=scaleyellowred,scale="none",key=TRUE, keysize=1.5,
            density.info="none", trace="none", cexRow=1,cexCol=1,margins = c(7, 7),lhei)
         )}
  


#' @title plots the 3D barplot
#' @description It shows a barplot of 5 networks given by the user with a summary representation of number of nodes, edges, and miRNAs (log values)
#' @param Edges_1net int number of edges in the 1 net
#' @param Edges_2net int number of edges in the 2 net
#' @param Edges_3net int number of edges in the 3 net
#' @param Edges_4net int number of edges in the 4 net
#' @param Edges_5net int number of edges in the 5 net
#' @param NODES_1net int number of nodes in the 1 net
#' @param NODES_2net int number of nodes in the 2 net
#' @param NODES_3net int number of nodes in the 3 net
#' @param NODES_4net int number of nodes in the 4 net
#' @param NODES_5net int number of nodes in the 5 net
#' @param nmiRNAs_1net int number of miRNAs in the 1 net
#' @param nmiRNAs_2net int number of miRNAs in the 2 net
#' @param nmiRNAs_3net int number of miRNAs in the 3 net
#' @param nmiRNAs_4net int number of miRNAs in the 4 net
#' @param nmiRNAs_5net int number of miRNAs in the 5 net
#' @examples 
#' SpidermiRvisualize_3Dbarplot(Edges_1net=1041003,Edges_2net=100016,Edges_3net=3008,
#' Edges_4net=1493,Edges_5net=1598,NODES_1net=16502,NODES_2net=13338,NODES_3net=1429,NODES_4net=675,
#' NODES_5net=712,nmiRNAs_1net=0,nmiRNAs_2net=74,nmiRNAs_3net=0,nmiRNAs_4net=0,nmiRNAs_5net=37)
#' @importFrom lattice cloud 
#' @importFrom  latticeExtra panel.3dbars
#' @export
#' @return barplot
SpidermiRvisualize_3Dbarplot<-function(Edges_1net,Edges_2net,Edges_3net,Edges_4net,Edges_5net,NODES_1net,NODES_2net,NODES_3net,NODES_4net,NODES_5net,nmiRNAs_1net,nmiRNAs_2net,nmiRNAs_3net,nmiRNAs_4net,nmiRNAs_5net){
  y <- as.factor(c("1net","2net","3net","4net","5net","1net","2net","3net","4net","5net","1net","2net","3net","4net","5net"))
  x<-c(round(log(Edges_1net)),round(log(Edges_2net)),round(log(Edges_3net)),round(log(Edges_4net)),round(log(Edges_5net)),round(log(NODES_1net)),round(log(NODES_2net)),round(log(NODES_3net)),round(log(NODES_4net)),round(log(NODES_5net)),round(log(nmiRNAs_1net)),round(log(nmiRNAs_2net)),round(log(nmiRNAs_3net)),round(log(nmiRNAs_4net)),round(log(nmiRNAs_5net)))
  z<-as.factor(c("edges","edges","edges","edges","edges","nodes","nodes","nodes","nodes","nodes","miRNAs","miRNAs","miRNAs","miRNAs","miRNAs"))
  cloud(x
        ~y+z, panel.3d.cloud=panel.3dbars, col.facet='green',  screen = list(z = -50, x = -30),
        xbase=0.8, ybase=0.3, scales=list(arrows=FALSE, col=1), 
        par.settings = list(axis.line = list(col = "transparent")))
}

#' @title Visualize results obtained by SpidermiR analysis with the direction of the interaction (pharmaco-gene and miRNA-gene)
#' @description Visualize the network
#' @param data The input data is a dataframe containing network data.
#' @examples 
#' miRNA_cNET <-data.frame(gA=c('hsa-let-7a','hsa-miR-141'),gB=c('FOXM1','CDK'),stringsAsFactors=FALSE)
#' SpidermiRvisualize_direction(data=miRNA_cNET)
#' @importFrom visNetwork visNetwork visInteraction visLegend %>%
#' @export
#' @return 3D graphic
SpidermiRvisualize_direction<-function(data){
  colnames(data) <- c("V1", "V2") 
  from=c(data$V1)
  to=c(data$V2)
  label=as.character()
  for (i in 1:nrow(data)){
    label[i]<-as.character(paste("Edge",i))
  }
  length=rep(100,nrow(data))
  edges2<-data.frame(from,to,label,length)
  edges2$arrows<-""
  edges2$arrows[(grep("hsa",edges2$from))]<-"to"
  edges2$arrows[(grep("[a-z]",edges2$to))]<- "from" #pharmaco
  edges2$dashes<-""
  edges2$dashes=rep("FALSE",nrow(data))
  edges2$titles<-""
  #titles=as.character()
  for (i in 1:nrow(data)){
    edges2$titles[i]<-as.character(paste("Edge",i))
  }
  edges2$smooth<-""
  edges2$smooth<-c(rep("FALSE",nrow(data)))
  edges2$shadow<-""
  edges2$shadow<-c(rep("FALSE",nrow(data)))
  #edges2 = data.frame(from, to,label,length,arrows,dashes,titles,smooth,shadow) 
  m<-c(data$V1)
  m2<-c(data$V2)
  s<-c(m,m2)
  id<- c(unique(s))
  label<-id
  
  #label=as.character()
  
  #for (i in 1:length(id)){
  # label[i]<-as.character(paste("Node",i))
  #}
  nodes2=data.frame(id,label)
  nodes2$group<-""
  
  nodes2$group<- replace(nodes2$group, nodes2$group == "", "gene")#gene 
  nodes2$group[(grep("[a-z]",nodes2$id))]<- "pharmaco" #pharmaco
  nodes2$group[grep("hsa" ,nodes2$id)]<- "mirna" #mirna
  nodes2$group[grep("orf" ,nodes2$id)]<- "gene" #gene 
  nodes2$group[grep("Spry2" ,nodes2$id)]<- "gene" 
  #group=as.character()
  #for (i in 1:length(id)){
  # group[i]<-as.character(paste("group",i))
  #}
  
  nodes2$value<-""
  value<-1:length(id)
  
  #for (i in 1:length(id)){
    #nodes2$value[i]<-as.character(paste("Edge",i))
    nodes2$value<-value
    
   # }
  
  nodes2$shape<-""
  nodes2$shape<- replace(nodes2$shape, nodes2$shape == "", "circle")#gene 
  nodes2$shape[(grep("[a-z]",nodes2$id))]<- "ellipse" #pharmaco
  nodes2$shape[grep("hsa" ,nodes2$id)]<- "box" #mirna
  nodes2$shape[grep("orf" ,nodes2$id)]<- "circle" #gene 
  nodes2$shape<-as.factor(nodes2$shape)
  nodes2$title<-""
  for (i in 1:length(id)){
    nodes2$title[i]<-as.factor(paste("<p><b>",i,"</b><br>Node !</p>"))
  }
  nodes2$color<-""
  nodes2$color<- replace(nodes2$color, nodes2$color == "", "lightblue")#gene 
  nodes2$color[(grep("[a-z]",nodes2$id))]<- "green" #pharmaco
  nodes2$color[grep("hsa" ,nodes2$id)]<- "orange" #mirna
  nodes2$color[grep("orf" ,nodes2$id)]<- "blue" #gene 
  nodes2$color<-as.factor(nodes2$color)
  
  nodes2$shadow<-""
  nodes2$shadow<-c(rep("FALSE",length(id)))
  
  nodes2$shadow<-as.logical(nodes2$shadow)
  
  
  return(visNetwork(nodes2, edges2[,c(1,2,5)], width = "200%")%>% visLegend(useGroups = FALSE, addNodes = data.frame(label = c("mirna","pharmaco","gene"), shape = c("box","ellipse","circle"),color=c("orange","green","lightblue") )
  )%>% visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) )
}



