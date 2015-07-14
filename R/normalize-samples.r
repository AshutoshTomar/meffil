#' Normalize Infinium HumanMethylation450 BeadChips
#'
#' Normalize a set of samples using their normalized quality control objects.
#'
#' @param norm.objects The list or sublist of \code{\link{meffil.normalize.quantiles}()}.
#' @param pseudo Value to add to the denominator to make the methylation
#' estimate more stable when calculating methylation levels (Default: 100).
#' @param just.beta If \code{TRUE}, then return just the normalized methylation levels; otherwise,
#' return the normalized methylated and unmethylated matrices (Default: TRUE).
#' @param cpglist.remove Optional list of CpGs to exclude from final output
#' @param verbose If \code{TRUE}, then detailed status messages are printed during execution (Default: \code{FALSE}).
#' @param ... Arguments passed to \code{\link[parallel]{mclapply}()}.
#' @return If \code{just.beta == TRUE}, the normalized matrix of 
#' methylation levels between between 0 and 1
#' equal to methylated signal/(methylated + unmethylated signal + pseudo).
#' Otherwise, a list containing two matrices, the normalized methylated and unmethylated signals.
#' 
#' @export
meffil.normalize.samples <- function(norm.objects, 
                                     pseudo=100,
                                     just.beta=T,
                                     cpglist.remove=NULL,
                                     filename=NULL,
                                     max.bytes=2^30-1, ## maximum number of bytes for mclapply
                                     verbose=F,
                                     ...) {
    stopifnot(length(norm.objects) >= 2)
    stopifnot(all(sapply(norm.objects, is.normalized.object)))

    sites <- meffil.get.sites()
    if(!is.null(cpglist.remove))
        sites <- setdiff(sites, cpglist.remove)
    
    ret <- mcsapply.safe(
        norm.objects,
        FUN=function(object) {
            ret <- meffil.normalize.sample(object, verbose=verbose)
            if (just.beta)
                get.beta(unname(ret$M[sites]), unname(ret$U[sites]), pseudo)
            else
                c(unname(ret$M[sites]), unname(ret$U[sites]))
        },
        ...,
        max.bytes=max.bytes)

    if (!just.beta) {
        ret <- list(M=ret[1:length(sites),],
                    U=ret[(length(sites)+1):nrow(ret),])
        rownames(ret$M) <- rownames(ret$U) <- sites
    }
    else
        rownames(ret) <- sites
    ret
}

#' We use \code{\link[parallel]{mclapply}()} to reduce running time by taking advantage of the fact
#' that each sample can be normalized independently of the others.
#' Unfortuantely \code{\link[parallel]{mclapply}()} has two important limitations.
#' The size in memory of the list returned may be at most 2Gb otherwise
#' \code{\link[parallel]{mclapply}()} fails with the following error message:
#'    Error in sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE)) :
#'    long vectors not supported ...
#' 
#' A non-elegant solution to this problem is to guess the size of each element
#' in the returned list and then apply \code{\link[parallel]{mclapply}()} sequentially to a sequence
#' appropriately sized input subsets.
#' A solution for \code{lapply} is to allocate the final object (e.g. a matrix)
#' prior to calling \code{lapply} and then populate the object during
#' the call to \code{lapply} using the global assignment operator '<<-'.
#' Unfortunately this is not a solution for \code{\link[parallel]{mclapply}()}
#' because \code{\link[parallel]{mclapply}()} immediately
#' duplicates the object, applies any modifications to the duplicate
#' and then deletes it prior to completion losing all modifications.
#' I'm not sure why the duplicate is not copied onto the original.
#' This is a solution if the object is a \code{\link{bigmemory::big.matrix}}.
#' However, we tried this but encountered random errors.
#' Sometimes the function completed without incident but other times,
#' with the same data, entire columns of the output matrix would be NA,
#' implying that the meffil.normalize.sample() function failed.
#' However, no errors were generated (tested with a tryCatch).
#' It seems that mclapply and big.matrix do not play well together all the time.
#' We have replaced this with the less elegant approach implemented in mcsapply().

