
#' mclapply
#'
#' mclapply without the 2Gb output memory bound.
#'
#' @param X Same as \code{\link[parallel]{mclapply}()}:
#' A vector (atomic or list) or an expressions vector.
#' Other objects (including classed objects) will be coerced by
#' \code{\link{as.list}()}.
#' @param FUN Same as \code{\link[parallel]{mclapply}()}:
#' The function to be applied to each element of 'X'.
#' @param ... Optional arguments to \code{FUN} and \code{\link[parallel]{mclapply}()}.
#' @param ret.bytes Number of bytes needed to store
#' each object returned by \code{FUN}, i.e. \code{FUN(X[[i]])}.
#' Number of bytes can be computed using the function \code{\link{object.size}()}.
#' If this value underestimates the actual memory requirements
#' and the returned list requires more than 2Gb of storage space,
#' then the function will fail with the following error message:
#'    Error in sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE)) :
#'    long vectors not supported ...
#' @param max.bytes The size in memory of the largest object that can
#' be returned by \code{\link[parallel]{mclapply}} (Default: 2Gb-1).
#' @return Same as \code{\link[parallel]{mclapply}()}:
#' a list of the same length as 'X' and named by 'X'.
#' Element i is equal to \code{FUN(X[[i]])}.
#'
#' @export
meffil.mclapply <- function (X, FUN, ..., ret.bytes=NA, max.bytes=2^30-1, tempdir=tempdir()) {
    stopifnot(!is.na(ret.bytes) & ret.bytes <= max.bytes)
    n.fun <- floor(max.bytes/ret.bytes)
    n.mclapply <- ceiling(length(X)/n.fun)
    partitions <- partition.integer.subsequence(1,length(X),n.mclapply)
    filenames <- sapply(1:nrow(partitions), function(i) tempfile(tmpdir=tempdir))
    for (i in 1:nrow(partitions)) {
        ret <- mclapply(X[partitions[i,"start"]:partitions[i,"end"]], FUN, ...)
        save(ret, file=filenames[i])
        rm(ret)
        gc()
    }
    do.call(c, lapply(filenames, function(filename) {
        msg("meffil.mclapply is loading", filename, "of", length(filenames), "filenames")
        load(filename, envir=sys.frame(sys.nframe()))
        unlink(filename)
        ret
    }))
}


partition.integer.subsequence <- function(start, end, n) {
    stopifnot(start <= end)
    stopifnot(n <= end-start+1)
    partitions <- floor(seq(start,end+1,length.out=n+1))
    cbind(start=head(partitions, n=-1),
          end=tail(partitions, n=-1) - 1)
}
