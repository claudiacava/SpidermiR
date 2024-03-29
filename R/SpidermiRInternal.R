#' @import utils
#' @importFrom httr GET content stop_for_status
.DownloadURL <-
  function(site){
    response <- GET(site)
    stop_for_status(response)
    x<-content(response,"text")
    x <- unlist(strsplit(x,"\n"))
    xi <- x[grep("href", x)]
    x <- sapply(strsplit(xi, "href"), function(y) y[2])
    x <- sapply(strsplit(x, ">"), function(y) y[2])
    return(x)
  }

#identifier_mappings for a particular 
.identifier<-function(organism_1) {
  Tableorganism<-paste(.url_cache$get("geneMania"),organism_1,"/identifier_mappings.txt",sep="")
  Tableorganism<-unlist(Tableorganism)
  return(Tableorganism)
}






#it creates link 
.url_cache <- local({
  env <- new.env(parent=emptyenv())
    #env[["miRwalk"]] <-
   # "http://zmf.umm.uni-heidelberg.de/apps/zmf/mirwalk2/downloads/vtm/hsa-vtm-gene.rdata.zip"
  env[["hmdd"]] <-
    "http://www.cuilab.cn/static/hmdd3/data/alldata.txt"
  env[["mirandola"]] <-
    "http://mirandola.iit.cnr.it/download/miRandola_version_02_2017.txt"
  env[["geneMania"]] <-
    "http://genemania.org/data/current/"
	  env[["dgidb"]] <-
	"https://www.dgidb.org/data/monthly_tsvs/2020-Sep/interactions.tsv"
	 # env[["matador"]] <-	  
	 # "http://matador.embl.de/media/download/matador.tsv.gz"
	 # env[["miRTarBase"]] <-
	  #  "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/cache/download/8.0/miRTarBase_SE_WR.xls"
	  env[["Tiger"]] <-
	    "http://bioinfo.wilmer.jhu.edu/tiger/download/hs2tissue-Table.txt"
	  env[["Tiger_mapping"]] <-
	    "http://bioinfo.wilmer.jhu.edu/tiger/download/symbol2hs-Table.txt"
	  env[["LNC"]] <-
	    "http://genome.igib.res.in/lncRNome/downloads/general_information.txt"
	  env[["LNC_mirna"]] <-
	  "http://genome.igib.res.in/lncRNome/downloads/miRNA_binding_sites.txt"
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

#internal function 
int<-function(at){
  se2=list()
  for (j in 1:nrow(at)){
    az<-at[j,]
    de<-as.data.frame(az)
    de[,2]<-rep(rownames(at)[j],length(de))
    de<-de[complete.cases(de),]
    se2[[j]]<-de

  }
  ds<-do.call("rbind", se2)
  ds<-ds[c(2,1)]
  # merging miRtar and miRNA walk information
  colnames(ds) <- c("V1", "V2")

  return(ds)
}

