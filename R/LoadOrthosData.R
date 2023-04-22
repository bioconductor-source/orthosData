#' Cache orthosData models
#'
#' Download in cache a set of orthosData keras models from ExperimentHub.
#'
#' @details
#' The function pre-caches a set of pre-trained orthosData keras models for a
#' given organism.
#' These are the models required to perform inference in `orthos`.
#' For each organism they are of three types:
#'
#' \strong{ContextEncoder}: The encoder component of a Variational Autoencoder
#' (VAE). Used to  produce a latent encoding of a given gene expression profile
#' (i.e context).
#'
#' - Input is a gene expression vector (shape=N, where N is the number or
#' `orthos` gene features) in the form of log2-transformed library normalized
#' counts (log2 counts per million, log2CPMs).
#'
#' - Output is a 64-d latent representation of the context.
#'
#' \strong{DeltaEncoder}: The encoder component of a conditional Variational
#' Autoencoder (cVAE).
#'
#' Used to produce a latent encoding of a contrast between two conditions (i.e
#' delta).
#'
#' - Input is a vector of gene expression contrasts (shape=N) in the form of
#' gene log2 CPM ratios (log2 fold changes, log2FCs),
#' concatenated with the corresponding context encoding.
#'
#' - Output is a 512-d latent representation of the contrast, conditioned on
#' the context.
#'
#' \strong{DeltaDecoder}: The decoder component of the same cVAE as above.
#' Used to produce the decoded version of the contrast between two conditions.
#'
#' - Input is the concatenated vector of the delta and context latent encodings.
#'
#' - Output is the decoded contrast vector (shape=N), conditioned on the context.
#'
#' For more details on model architecture and use of these models in `orthos`
#' please refer to the `orthos` package vignette:
#' vignette("orthosIntro", package = "orthos").
#'
#'
#' @param organism Character scalar selecting the organism for which to load the
#'     contrast database. One of \code{"Human"} or \code{"Mouse"}.
#' @param ARCH4v Version of ARCHS4 used to build the contrastDB.
#' @param verbose Logical scalar indicating whether to print messages along
#'     the way.
#'
#' @export
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom AnnotationHub query hubUrl getInfoOnIds
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom stringr str_extract
#'
#' @examples
#' \dontrun{
#' GetorthosModels(organism = "Mouse")
#' }
GetorthosModels <- function(organism = c("Human", "Mouse"),
                            ARCHS4v = "v212",
                            verbose = TRUE)
{
    organism <- tolower(match.arg(organism))
    hub <- ExperimentHub::ExperimentHub()

    ## Models info from ExperimentHub:
    query_keys <- c("orthosData", "coder", organism, "ARCHS4", ARCHS4v)
    hub_query_results <- AnnotationHub::query(hub, query_keys)

    if (verbose) {
        message(paste("caching resources: ", hub_query_results$title,
                      collapse = "\n"))
        AnnotationHub::cache(hub_query_results)
    } else {
        suppressMessages({
            AnnotationHub::cache(hub_query_results)
        })
    }
}



#' Cache a web resource to a BiocFileCache location
#'
#' Cache a web resource to a BiocFileCache location
#' typically at: \code{ExperimentHub::getExperimentHubOption("CACHE")}.
#'
#' @param url url link of the resource file to be cached.
#' @param rname resource name
#' @param fname Options are ‘unique’ or ‘exact’. See  \code{BiocFileCache}
#'
#' @return the \code{dirname} of the cached objects.
#'
#' @author Panagiotis Papasaikas
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfccount bfcadd
#'     bfcneedsupdate bfcdownload
#' @importFrom ExperimentHub getExperimentHubOption

#'
#' @keywords internal
#' @noRd
.addToCache <- function(cache_path = ExperimentHub::getExperimentHubOption("CACHE"),
                        url = "",
                        rname = "",
                        fname = c("exact", "unique")) {

    fname <- match.arg(fname)
    bfc <- BiocFileCache::BiocFileCache(cache_path)

    # check if url is being tracked
    res <- BiocFileCache::bfcquery(bfc, url)

    if (BiocFileCache::bfccount(res) == 0L) {
        # if not in cache, add it
        message(paste("Caching: ", rname, collapse = "\n"))
        ans <- BiocFileCache::bfcadd(bfc, rname = rname, fpath = url,
                                     fname = fname)
        message(paste("File cached at: ", ans, collapse = "\n"))
    } else {
        # if it is in cache, get path to load
        ans <- res$rpath[1]
        message(paste(rname," already present in cache at:", ans,
                      collapse = "\n"))

        # check to see if the resource needs to be updated
        check <- BiocFileCache::bfcneedsupdate(bfc, res$rid[1])
        # check can be NA if it cannot be determined, choose how to handle
        if (is.na(check)) check <- TRUE
        if (check) {
                message(paste("updating resource", ans, collapse = "\n"))
            ans <- BiocFileCache::bfcdownload(bfc, res$rid[1])
        }
    }
# Return caching directory
ans <- dirname(ans)
return(ans)
}


#' Cache an orthosData contrast DB
#'
#' Download in cache HDF5 and RDS component files for an othosData contrast
#' database from ExperimentHub.
#' As these are components of a single HDF5SummarizedExperiment object they
#' HAVE to be cached with the exact prefix used at creation time (see
#' HDF5Array::saveHDF5SummarizedExperiment).
#'
#' @details The orthosData  contrast database contains over 100,000
#' differential gene expression experiments
#' compiled from the ARCHS4 database* of publicly available expression data.
#' Each entry in the database corresponds to a pair of RNAseq samples contrasting
#' a treatment vs a control condition.
#'
#' A combination of metadata-semantic and quantitative analyses was used to
#' determine the proper assignment of samples to such pairs in `orthosData`.
#'
#' The ~20,000 gene features/organism used in the database are "sanctioned"
#' according to several criteria (located on canonical chromosomes, no
#' pseudogenes, no ribosomal protein genes, detected in at least a small
#' fraction of the ARCHS4 database).
#'
#' The orthosData  contrast database contains assays with the original
#' contrasts in the form of gene expression log2 CPM ratios (i.e log2 fold
#' changes, log2FCs), precalculated, decoded and residual components of those
#' contrasts using the orthosData models as well as the gene expression context
#' of those contrasts in the form of log2-transformed library normalized counts
#' (i.e log2 counts per million, log2CPMs).
#' It also contains extensive annotation on both the `orthos` feature genes
#' and the contrasted conditions.
#'
#' For each organism the DB is stored as an HDF5SummarizedExperiment with an
#' HDF5 component that contains the gene assays and an rds component that
#' contains gene annotation in the rowData and the contrast annotation in the
#' colData.
#'
#' Note that because of the way that HDF5 datasets and serialized
#' SummarizedExperiments are linked in an HDF5SummarizedExperiment, the two
#' components -although relocatable- need to have the exact same filenames
#' as those used at creation time. In other words the files can be moved (or
#' copied) to a different directory or to a different machine and they will
#' retain functionality as long as both live in the same directory and are
#' never renamed.
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
#' @return the \code{dirname} of the cached objects.
#' @export
#' @author Panagiotis Papasaikas
#'
#' @importFrom AnnotationHub query hubUrl getInfoOnIds
#' @importFrom ExperimentHub ExperimentHub getExperimentHubOption
#' @importFrom stringr str_extract
#'
#' @examples
#' \dontrun{
#' GetorthosContrastDB(organism = "Mouse", mode="DEMO")
#'
#' se <- HDF5Array::loadHDF5SummarizedExperiment(dir = ExperimentHub::getExperimentHubOption("CACHE"),
#' prefix = "mouse_v212_NDF_c100_DEMO")
#'
#' }
#' @references{
#' *Lachmann, Alexander, et al.
#' "Massive mining of publicly available RNA-seq data from human and mouse."
#'  Nature communications 9.1 (2018): 1366
#' }
#'
GetorthosContrastDB <- function(organism = c("Human", "Mouse"),
                         mode = c("ANALYSIS", "DEMO"),
                         ARCHS4v = "v212",
                         verbose = TRUE)
{
    organism <- tolower(match.arg(organism))
    DEMO <- ifelse(mode == "DEMO", "_DEMO", "")
    hub <- ExperimentHub()

    ## RDS info from ExperimentHub:
    rds_rdatapath <- paste0(organism, "_", ARCHS4v, "_NDF_c100", DEMO, "se.rds")
    rds_hubObj <- AnnotationHub::query(hub, c("orthosData", rds_rdatapath))
    hubUrl <- AnnotationHub::hubUrl(rds_hubObj)

    info <- AnnotationHub::getInfoOnIds(rds_hubObj)
    title <- info$title
    fetch_id <- info$fetch_id
    fetch_url <- paste0(hubUrl, "/fetch/", fetch_id)
    fetch_url_real <- stringr::str_extract(grep(curlGetHeaders(fetch_url),
                                                pattern = "Location",
                                                value = TRUE),
                                           pattern = "http.*") # resolve redirects
    if (verbose) {
        .addToCache(url = fetch_url_real, rname = title)
    } else {
        suppressMessages({
            .addToCache(url = fetch_url_real, rname = title)
        })
    }


    ## HDF5 info from ExperimentHub:
    h5_rdatapath <- paste0(organism, "_", ARCHS4v, "_NDF_c100", DEMO,
                           "assays.h5")
    h5_hubObj <- query(hub, c("orthosData", h5_rdatapath))

    info <- AnnotationHub:::getInfoOnIds(h5_hubObj)
    title <- info$title
    fetch_id <- info$fetch_id
    fetch_url <- paste0(hubUrl, "/fetch/", fetch_id)
    fetch_url_real <- stringr::str_extract(grep(curlGetHeaders(fetch_url),
                                                pattern = "Location",
                                                value = TRUE),
                                           pattern = "http.*") # resolve redirects
    if (verbose) {
        cache_dir <- .addToCache(url = fetch_url_real, rname = title)
    } else {
        suppressMessages({
            cache_dir <- .addToCache(url = fetch_url_real, rname = title)
        })
    }
    cache_dir
}
