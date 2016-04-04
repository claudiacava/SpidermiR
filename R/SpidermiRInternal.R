

#' @import utils
#' @importFrom httr GET content stop_for_status
.DownloadURL <-
  function(Site_g){
    response <- GET(Site_g)
    stop_for_status(response)
    x<-content(response,"text")
    x <- unlist(strsplit(x,"\n"))
    xi <- x[grep("href", x)]
    x = sapply(strsplit(xi, ">"), function(y) y[2])
    return(x)
  }

#identifier_mappings for a particular 
.identifier<-function(organism_1) {
  Tableorganism<-paste(.url_cache$get("geneMania"),organism_1,"/identifier_mappings.txt",sep="")
  Tableorganism<-unlist(Tableorganism)
  return(Tableorganism)
}




#' @importFrom networkD3 simpleNetwork
.SpidermiRvisualize_gene<-function(data){
  s<-as.data.frame(data)
  Source <- (s[1])
  Target <- (s[2])
  NetworkData <- data.frame(Source, Target)
  simpleNetwork(NetworkData,linkColour = "gray",textColour = "black",zoom = TRUE)
}


#it creates link 
.url_cache <- local({
  env <- new.env(parent=emptyenv())
  env[["miRtar"]] <-
    "http://watson.compbio.iupui.edu:8080/miR2Disease/download/miRtar.txt"
  env[["miRwalk"]] <-
    "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/hsa-vtm-gene.rdata.zip"
  env[["miR2Disease"]] <-
    "http://watson.compbio.iupui.edu:8080/miR2Disease/download/AllEntries.txt"
  env[["mirandola"]] <-
    "http://atlas.dmi.unict.it/mirandola/download/miRandola_version_1.7_10-06-2015.csv"
  env[["geneMania"]] <-
    "http://genemania.org/data/current/"
  list(
    get=function(elt) {
      stopifnot(is.character(elt), length(elt) == 1L, elt %in%
                  ls(env))
      env[[elt]]
    }, set=function(elt, url) {
      stopifnot(is.character(elt), length(elt) == 1L, elt %in%
                  ls(env))
      stopifnot(is.character(url), length(url) == 1L)
      env[[elt]] <- url
      env[[elt]]
    })
})




