#' @title Download the network from GeneMania. 
#' @description
#'      SpidermiRdownload_net function will download the data 
#'
#' @param data The SpidermiRquery_spec_networks output
#' @examples
#' org<-SpidermiRquery_species(species)
#' net_shar_prot<-SpidermiRquery_spec_networks(organismID = org[18,],
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

