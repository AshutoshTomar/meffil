#' Save detection p-value matrix to GDS file
#'
#' @param qc.objects A list of outputs from \code{\link{meffil.create.qc.object}()}.
#' @param gds.filename If not \code{NULL} (default), then saves the output to a GDS (Genomic Data Structure).
#' This is for cases where the output is too large to fit into main memory.
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return Matrix of probe detection p-values.
#' If \code{gds.filename} is not \code{NULL}, then
#' the output is saved to the GDS file rather than retained
#' in memory and returned.
#' The library 'gdsfmt' must be installed in this case.
#' 
#' @export
meffil.save.detection.pvalues <- function(qc.objects,
                                          gds.filename=NULL,
                                          max.bytes=2^30-1,
                                          verbose=F,
                                          ...) {
    stopifnot(all(sapply(qc.objects, is.qc.object)))

    featuresets <- sapply(qc.objects, function(qc.object) qc.object$featureset)
    featureset <- featuresets[1]

    if (is.list(featuresets)) ## backwards compatibility
        featureset <- featuresets <- "450k"

    if (any(featuresets != featureset)) 
        stop("Multiple feature sets were used to create these QC objects:",
             paste(unique(featuresets), collapse=", "))

    feature.names <- meffil.get.features(featureset)$name
    
    if (!all(sapply(qc.objects, function(qc.object) exists.rg(qc.object$basename))))
         stop("IDAT files are not accessible for all QC objects")

    require(gdsfmt)
    mcsapply.to.gds(
        qc.objects,
        function(object) {
            if (is.null(object$featureset)) ## backwards compatibility
                object$chip <- "450k" 
            rg <- read.rg(object$basename, verbose=verbose)        
            probes <- meffil.probe.info(object$chip)        
            pvalues <- extract.detection.pvalues(rg, probes, verbose=verbose)
            pvalues[feature.names]        
        },
        ...,
        gds.filename=gds.filename,
        storage="float64",
        max.bytes=max.bytes)
    gds.filename
}


#' Retrieve detection p-values from GDS file
#'
#' @param gds.filename Name of GDS file generated by \code{\link{meffil.save.detection.pvalues}()}.
#' @param sites Names of CpG sites to load, if `NULL` then load all (Default: NULL).
#' @param samples Names of samples to load, if `NULL` then load all (Default: NULL).
#' @return Matrix of methylation levels with rows corresponding to CpG sites
#' and columns to samples.  Rows restricted \code{sites} if not \code{NULL},
#' and columns restricted to \code{samples} if not \code{NULL}.
#' 
#'
#' @export
meffil.gds.detection.pvalues <- function(gds.filename, sites=NULL, samples=NULL) {
    retrieve.gds.matrix(gds.filename, sites, samples)
}
