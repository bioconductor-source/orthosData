#' Cache a web resource to a BiocFileCache location 
#'
#' Cache a web resource to a BiocFileCache location 
#' typically at: \code{Sys.getenv("EXPERIMENT_HUB_CACHE")}.
#'
#' @param url url link of the resource file to be cached. 
#' @param rname resource name
#' @param fname Options are ‘unique’ or ‘exact’. See  \code{BiocFileCache}
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfccount bfcadd bfcneedsupdate bfcdownload
#'
#' @keywords internal
#' @noRd
.addToCache <- function(cache_path=Sys.getenv("EXPERIMENT_HUB_CACHE"),
                        url="",
                        rname="",
                        fname=c("unique", "exact")) {
    
    bfc <- BiocFileCache::BiocFileCache(cache_path)
    
    # check if url is being tracked
    res <- BiocFileCache::bfcquery(bfc, url)
    
    if (BiocFileCache::bfccount(res) == 0L) {
        # if not in cache, add it
        message(paste("Caching: ", rname, collapse = "\n"))
        ans <- BiocFileCache::bfcadd(bfc, rname=rname, fpath=url, fname = fname)
        message(paste("File cached at: ", ans, collapse = "\n"))
        
    } else {
        # if it is in cache, get path to load
        ans <- res$rpath[1]
        message(paste(rname," already present in cache at:", ans, collapse = "\n"))
        
        # check to see if the resource needs to be updated
        check <- BiocFileCache::bfcneedsupdate(bfc, res$rid[1])
        # check can be NA if it cannot be determined, choose how to handle
        if (is.na(check)) check <- TRUE
        if (check){
            if (verbose) {
                message(paste("updating resource", ans, collapse = "\n"))
            }
            ans <- BiocFileCache::bfcdownload(bfc, res$rid[1])
        }
    }
    
}


#' Cache an orthosData contrast DB
#'
#' Download in cache HDF5 and RDS component files for an othosData contrast DB from ExperimentHub.
#' As these are components of a single HDF5SummarizedExperiment object they HAVE to be cached
#' with the exact prefix used at creation time (see HDF5Array::saveHDF5SummarizedExperiment)
#'
#'
#' @param organism Character scalar selecting the organism for which to load the
#'     contrast database. One of \code{"Human"} or \code{"Mouse"}.
#' @param mode When in "ANALYSIS" mode (default) the complete contrast DB is
#'     cached. "DEMO" mode caches a small "toy" database for the queries.
#'     "DEMO" should only be used for testing/demonstration purposes
#'     and never for actual analysis purposes.
#' @param ARCH4v Version of ARCHS4 used to build the contrastDB.     
#' @param verbose Logical scalar indicating whether to print messages along
#'     the way.
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom AnnotationHub query hubUrl getInfoOnIds
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom stringr str_extract
#' 
#' @examples
#' \dontrun{
#' GetorthosContrastDB(organism = "Mouse", mode="DEMO")
#' }
GetorthosContrastDB <- function(organism = c("Human","Mouse"),
                         mode = c("ANALYSIS", "DEMO"),
                         ARCHS4v = "v212",
                         verbose=TRUE)
{
    organism <- tolower(match.arg(organism))
    DEMO <- ifelse(mode=="DEMO","_DEMO","")
    hub <- ExperimentHub()

        ## RDS info from ExperimentHub:
        rds_rdatapath <- paste0(organism,"_",ARCHS4v,"_NDF_c100",DEMO,"se.rds" )
        rds_hubObj <- AnnotationHub::query(hub, rds_rdatapath)
        hubUrl <- AnnotationHub::hubUrl(rds_hubObj )
        
        info <- AnnotationHub::getInfoOnIds(rds_hubObj )
        title <- info$title 
        fetch_id <- info$fetch_id
        fetch_url <- paste0( hubUrl,"/fetch/", fetch_id)
        fetch_url_real <- stringr::str_extract(grep(curlGetHeaders(fetch_url), pattern = "Location", value = T), pattern = "http.*") # resolve redirects
        if(verbose){
        .addToCache(url=fetch_url_real,rname=title)
        }
        else{
            suppressMessages({
                .addToCache(url=fetch_url_real,rname=title)    
                })
        }
            
        
        ## HDF5, from ExperimentHub:
        h5_rdatapath <- paste0(organism,"_",ARCHS4v,"_NDF_c100",DEMO,"assays.h5" )
        h5_hubObj <- query(hub, h5_rdatapath )
        
        info <- AnnotationHub:::getInfoOnIds(h5_hubObj )
        title <- info$title 
        fetch_id <- info$fetch_id
        fetch_url <- paste0( hubUrl,"/fetch/", fetch_id)
        fetch_url_real <- stringr::str_extract(grep(curlGetHeaders(fetch_url), pattern = "Location", value = T), pattern = "http.*") # resolve redirects
        if(verbose){
            .addToCache(url=fetch_url_real,rname=title)
        }
        else{
            suppressMessages({
                .addToCache(url=fetch_url_real,rname=title)    
            })
        }

}









#' Cache orthosData models
#'
#' Download in cache a set of orthosData keras models from ExperimentHub.
#'
#' @param organism Character scalar selecting the organism for which to load the
#'     contrast database. One of \code{"Human"} or \code{"Mouse"}.
#' @param ARCH4v Version of ARCHS4 used to build the contrastDB.     
#' @param verbose Logical scalar indicating whether to print messages along
#'     the way.
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom AnnotationHub query hubUrl getInfoOnIds
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom stringr str_extract
#' 
GetorthosModels <- function(organism = c("Human","Mouse"),
                                ARCHS4v = "v212",
                                verbose=TRUE)
{
    organism <- tolower(match.arg(organism))
    hub <- ExperimentHub::ExperimentHub()
    
    ## Model info from ExperimentHub:
    query_keys <- c( "coder",organism, "ARCHS4", ARCHS4v   )
    hub_query_results <- AnnotationHub::query(hub, query_keys)
    
    if(verbose){
        message(paste("caching resources: ",hub_query_results$title, collapse = "\n"))
        AnnotationHub::cache(hub_query_results)
    }
    else{
        suppressMessages({
            AnnotationHub::cache(hub_query_results)    
        })
    }
}


