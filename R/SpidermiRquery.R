#' @title Searching  by network species
#' @description The user can visualize the species supported by GeneMania, using the function SpidermiRquery_species .
#' @param species a variable parameter 
#' @examples org<-SpidermiRquery_species(species)
#' @export
#' @return table of species
SpidermiRquery_species <- function(species) {
  site <- .url_cache$get("geneMania")
  site1<-.DownloadURL(site)
  site1<-site1[-1]
  site1<-site1[-1]
  Tableorganism <- matrix(0, length(site1), 1)
  rownames(Tableorganism)<-paste("organism",c(1:length(site1)),sep="")
  colnames(Tableorganism) <- c("link")
  organism_2=list()
  for ( tbi in 1:nrow(Tableorganism)){
    organism <- site1[tbi]
    #qst_find_sitep <-gsub("<a href=","", as.matrix(unlist(strsplit(organism,">")))[1])
    qst_find_sitep <-gsub("</a","", as.matrix(unlist(strsplit(organism,">")))[1])
    #qst_find_site2_subp <- substr(qst_find_sitep,3,nchar(qst_find_sitep)-1)
    organism_2[tbi]<-qst_find_sitep
    newsite_tofind <- paste(site,qst_find_sitep,sep="")
    Tableorganism[tbi,"link"] <- newsite_tofind
  }
  organismID<-gsub("/","",organism_2)
  rownames(Tableorganism)<-organismID
  tabOrgd <- matrix(0,nrow(as.matrix(organismID)),2)
  colnames(tabOrgd) <- c("Number","Species")
  tabOrgd <- as.data.frame(tabOrgd)
  tabOrgd$Number <- rownames(organismID)
  tabOrgd$Species <- organismID
 # organismID  <- as.matrix(organismID)
  #colnames(organismID) <- "Species"
  tabOrgd<-tabOrgd[- grep("COMBINED", tabOrgd$Species),]
  tabOrgd<-as.data.frame(tabOrgd)
  tabOrgd<-tabOrgd[- grep("README", tabOrgd$tabOrgd),]
  tabOrgd<-as.data.frame(tabOrgd)
  return(tabOrgd)
  
}






#' @title Network categories
#' @description The user can visualize the network types supported by GeneMania for a specific specie using SpidermiRquery_networks_type
#' @param organismID describes index of a specific specie obtained by SpidermiRquery_species output
#' @export
#' @examples
#' org<-SpidermiRquery_species(species)
#' net_type<-SpidermiRquery_networks_type(organismID=org[9,])
#' @return a list of network categories in a specie indicated.
SpidermiRquery_networks_type<-function(organismID) {
  vd<-list()
  Tableorganism<-paste(.url_cache$get("geneMania"),organismID,"/networks.txt",sep="")
  sd<-read.delim(Tableorganism,header = TRUE,stringsAsFactors=FALSE)
  net_t<-sd[!duplicated(sd$Network_Group_Name), ]
  net_t<-net_t$Network_Group_Name
  return(net_t)
}











#' @title Searching by network categories
#' @description The user can visualize the database or reference where the information came from
#' @param organismID describes index of a specific specie obtained by SpidermiRquery_species output
#' @param network The network type the user is interested in.  
#' Example:
#' \tabular{ll}{
#'COexp \tab   Co-expression \cr
#'PHint  \tab   Physical_interactions \cr
#'COloc  \tab   Co-localization \cr
#'GENint  \tab   Genetic_interactions \cr
#'PATH  \tab   Pathway \cr
#'SHpd  \tab   Shared_protein_domains \cr
#'}
#' @export
#' @examples
#' org<-SpidermiRquery_species(species)
#' net_shar_prot<-SpidermiRquery_spec_networks(organismID = org[9,],
#'                                         network = "SHpd")
#'                                         
#'                                         
#' @return a list of the database or reference where the information came from.
SpidermiRquery_spec_networks<-function(organismID,network) {
  vd<-list()
  Tableorganism<-paste(.url_cache$get("geneMania"),organismID,"/",sep="")
  Site_sec <- .DownloadURL(Tableorganism)
  #qst_find_sitep_sub<-gsub("<a href=","",Site_sec)
  qst_find_sitep <-gsub("</a","", as.matrix(unlist(strsplit(Site_sec,">"))))
  #qst_find_site2_subp_sub <- substr(qst_find_sitep_sub,3,nchar(qst_find_sitep_sub)-1)
  #qst_find_site2_subp_sub<-qst_find_site2_subp_sub[-1]
  newsite_tofind_sub_sub <- list(paste(Tableorganism,qst_find_sitep,sep=""))
  mylist<-newsite_tofind_sub_sub
  for (i in 1:length(mylist)){
    mylist_gm<-mylist[i]
    x <- unlist(mylist_gm)
    #Co-expression
    if(network == "COexp"){
      x<-x[grep("/Co-expression", x)]}
    #Physical_interactions
    if(network == "PHint"){
      x<-x[grep("/Physical_interactions", x)]}
    #Co-localization
    if(network == "COloc"){
      x<-x[grep("/Co-localization", x)]}
    #Genetic_interactions
    if(network == "GENint"){
      x<-x[grep("/Genetic_interactions", x)]}
    #pathway
    if(network == "PATH"){
      x<-x[grep("/Pathway", x)]}
    #Shared_protein_domains
    if(network == "SHpd"){
      x<-x[grep("/Shared_protein_domains", x)]}
    vd[[i]]<-x
    vd<-unlist(vd)
  }
  return (vd) 
}





#' @title Visualize disease categories
#' @description The user can visualize the disease supported by SpidermiR
#' @param diseaseID variable name
#' @export
#' @examples
#' disease<-SpidermiRquery_disease(diseaseID)
#' @return a list of disease.
SpidermiRquery_disease<-function(diseaseID) {
all_entries<-.url_cache$get("miR2Disease")
disease_ref<-read.delim(all_entries,header = FALSE,quote = "",stringsAsFactors=FALSE)
disease<- unique(disease_ref$V2)
disease<-disease[disease != ""]
return (disease)
}




