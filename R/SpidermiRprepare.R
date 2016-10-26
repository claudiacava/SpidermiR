#' @title Prepare matrix of gene network from Genamania with Ensembl Gene ID, and gene symbols
#' @description The user in this step obtained a gene network matrix with the integration of gene symbols ID.
#' @return A list of tables.
#' @param organismID is the index of SpidermiRquery_spec_networks output
#' @param data is the output of function SpidermiRdownload_net
#' @examples
#' geneSymb_net<-SpidermiRprepare_NET(organismID = org[9,],
#'                                        data = out_net)
#' @export


SpidermiRprepare_NET <- function(organismID,data){
    #prova<-data

  #find an identifier_mappings for a species
  organismID <-.identifier(organismID)
  
    filename_id<-read.delim(organismID,header = TRUE,stringsAsFactors=FALSE)
    filename_id_gs<-filename_id[filename_id$Source=="Gene Name",]
    list_r<-list()
    for(i in 1:length(data)){
      print(paste("Preprocessing of the network n. ", i, " of ", length(data), sep = ""))     
      a<-data[[i]]
      am<-match(a$Gene_A,filename_id_gs$Preferred_Name)
      st3 <- cbind( a, filename_id_gs[am,] )
      st4<-na.omit(st3)
      colnames(st4)[5]<-c("gene_symbolA")
      sed<-match(st4$Gene_B,filename_id_gs$Preferred_Name)
      st5 <- cbind( st4, filename_id_gs[sed,] )
      st6<-na.omit(st5)
      colnames(st6)[8]<-c("gene_symbolB")
      list_r[[i]]<-st6

    }
    return(list_r)
  }