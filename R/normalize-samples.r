#' Normalize Infinium HumanMethylation450 BeadChips
#'
#' Normalize a set of samples using their normalization objects.
#'
#' @param objects A list or sublist returned by \code{\link{meffil.normalize.objects}()}.
#' @param beta If \code{TRUE} (default), then the function returns
#' the normalized matrix of methylation levels; otherwise, it returns
#' the normalized matrices of methylated and unmethylated signals.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}
#' except for \code{ret.bytes}.
#' @return Matrix of normalized methylation levels if \code{beta} is \code{TRUE};
#' otherwise matrices of normalized methylated and unmethylated signals.
#' Matrices returned have one column per sample and one row per CpG site.
#' Methylation levels are values between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#'
#' @export
meffil.normalize.samples <- function(objects, beta=T, pseudo=100,
                                     probes=meffil.probe.info(), ...) {
    stopifnot(length(objects) >= 2)

    n.sites <- length(unique(probes$name))
    if (beta) {
        ret.bytes <- object.size(rep(NA_real_, n.sites))
        do.call(cbind, meffil.mclapply(objects, function(object) {
            meffil.get.beta(meffil.normalize.sample(object, probes=probes), pseudo)
        }, ret.bytes=ret.bytes, ...))
    }
    else {
        ret.bytes <- object.size(list(rep(NA_real_, n.sites),
                                      rep(NA_real_, n.sites)))
        ## approx 8 x n.sites x 2
        ret <- meffil.mclapply(objects,
                               meffil.normalize.sample, probes=probes,
                               ret.bytes=ret.bytes, ...)
        list(M=sapply(ret, function(x) x$M),
             U=sapply(ret, function(x) x$U))
    }
}
