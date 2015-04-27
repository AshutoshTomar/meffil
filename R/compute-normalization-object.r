#' Normalization object
#'
#' Create a normalization object for a given Infinium HumanMethylation450 BeadChip.
#'
#' @param samplesheet Data frame containing IDAT file and sample info (see \code{\link{meffil.read.samplesheet}} pr \code{\link{meffil.create.samplesheet}}).
#' @param number.quantiles Number of quantiles to compute for probe subset (Default: 500).
#' @param dye.intensity Reference intensity for scaling each color channel (Default: 5000).
#' @param probes Probe annotation used to construct the control matrix
#' (Default: \code{\link{meffil.probe.info}()}).
#' @param verbose If \code{TRUE}, then status messages are printed during execution (Default: \code{FALSE}).
#' @param detection.threshold Default value = 0.01.
#' All probes above this detection threshold detected.
#' @param bead.threshold Default value = 3.
#' All probes with less than this number of beads detected.
#' @return List containing control probe information, probe summaries
#' and quantiles.
#'
#' @export
meffil.compute.normalization.object <- function(samplesheet.row,
                                                number.quantiles=500,
                                                dye.intensity=5000,
                                                probes=meffil.probe.info(),
                                                verbose=F,
                                                detection.threshold=0.01,
                                                bead.threshold=3,
                                                sex.cutoff=-2) {
    stopifnot(number.quantiles >= 100)
    stopifnot(dye.intensity >= 100)
    stopifnot(nrow(probes) > 100000)
    stopifnot(samplesheet.row$Sex %in% c(NA, "F", "M"))

    rg <- meffil.read.rg(basename, probes, verbose=verbose)

    bad.probes.detectionp <- identify.bad.probes.detectionp(rg, detection.threshold, probes, verbose=verbose)

    bad.probes.beadnum <- identify.bad.probes.beadnum(rg, bead.threshold, probes, verbose=verbose)

    snp.probes <- get.snp.probes(rg, probes, verbose=verbose)

    controls <- extract.controls(rg, probes, verbose=verbose)

    rg <- meffil.background.correct(rg, probes, verbose=verbose)

    intensity.R <- calculate.intensity.R(rg, probes)
    intensity.G <- calculate.intensity.G(rg, probes)

    rg <- meffil.dye.bias.correct(rg, dye.intensity, probes, verbose=verbose)

    mu <- meffil.rg.to.mu(rg, probes, verbose=verbose)

    probes.x <- probes$name[which(probes$chr == "chrX")]
    x.signal <- median(log(mu$M[probes.x] + mu$U[probes.x], 2), na.rm=T)

    probes.y <- probes$name[which(probes$chr == "chrY")]
    y.signal <- median(log(mu$M[probes.y] + mu$U[probes.y], 2), na.rm=T)

    probs <- seq(0,1,length.out=number.quantiles)

    quantiles <- lapply(get.quantile.probe.subsets(probes), function(sets) {
        list(M=unname(quantile(mu$M[sets$M], probs=probs,na.rm=T)),
             U=unname(quantile(mu$U[sets$U], probs=probs,na.rm=T)))
    })

    msg("predicting sex", verbose=verbose)
    xy.diff <- y.signal-x.signal
    predicted.sex <- ifelse(xy.diff < sex.cutoff, "F","M")


    list(origin="meffil.compute.normalization.object",
         basename=basename,
         controls=controls,
         quantiles=quantiles,
         dye.intensity=dye.intensity,
         intensity.R=intensity.R,
         intensity.G=intensity.G,
         x.signal=x.signal,
         y.signal=y.signal,
         xy.diff=xy.diff,
         sex=samplesheet$Sex,
         sex.cutoff=sex.cutoff,
         predicted.sex=predicted.sex,
         median.m.signal=median(mu$M,na.rm=T),
         median.u.signal=median(mu$U,na.rm=T),
         bad.probes.detectionp=bad.probes.detectionp,
         bad.probes.beadnum=bad.probes.beadnum,
         bad.probes.detectionp.threshold=detection.threshold,
         bad.probes.beadnum.threshold=bead.threshold,
         snp.probes=snp.probes



         object$sex.cutoff <- sex.cutoff
         object$xy.diff <- xy.diff[i]
         object$predicted.sex <- predicted.sex[i]
         object$sex <- sex[i]


         )
}

is.normalization.object <- function(object) {
    (all(c("quantiles","dye.intensity","origin","basename","x.signal","y.signal","controls",
           "intensity.R","intensity.G")
         %in% names(object))
     && object$origin == "meffil.compute.normalization.object")
}

get.quantile.probe.subsets <- function(probes=meffil.probe.info()) {
    rm.na <- function(x) {
        x[which(is.na(x))] <- F
        x
    }

    is.iG <- rm.na(probes$type == "i" & probes$dye == "G")
    is.iR <- rm.na(probes$type == "i" & probes$dye == "R")
    is.ii <- rm.na(probes$type == "ii")
    is.genomic <- !is.na(probes$chr)
    is.sex <- rm.na(is.genomic & probes$chr %in% c("chrX","chrY"))
    is.x <- rm.na(is.genomic & probes$chr == "chrX")
    is.y <- rm.na(is.genomic & probes$chr == "chrY")
    is.autosomal <- rm.na(is.genomic & !is.sex)
    is.not.y <- rm.na(is.genomic & probes$chr != "chrY")

    get.probe.subsets <- function(in.subset) {
        list(M=probes$name[which(probes$target == "M" & in.subset)],
             U=probes$name[which(probes$target == "U" & in.subset)])
    }

    list(genomic.iG = get.probe.subsets(is.iG & is.genomic),
         genomic.iR = get.probe.subsets(is.iR & is.genomic),
         genomic.ii = get.probe.subsets(is.ii & is.genomic),
         autosomal.iG = get.probe.subsets(is.iG & is.autosomal),
         autosomal.iR = get.probe.subsets(is.iR & is.autosomal),
         autosomal.ii = get.probe.subsets(is.ii & is.autosomal),
         not.y.iG = get.probe.subsets(is.iG & is.not.y),
         not.y.iR = get.probe.subsets(is.iR & is.not.y),
         not.y.ii = get.probe.subsets(is.ii & is.not.y),
         sex = get.probe.subsets(is.sex),
         chry = get.probe.subsets(is.y),
         chrx = get.probe.subsets(is.x))
}

sex.specific.quantile.probe.subsets <- function() {
    c("genomic.iG",
      "genomic.iR",
      "genomic.ii",
      "not.y.iG",
      "not.y.iR",
      "not.y.ii",
      "sex",
      "chry",
      "chrx")
}

applicable.quantile.probe.subsets <- function(sex, both.sexes) {
    if (both.sexes && sex == "M") return(c("autosomal.iG","autosomal.iR","autosomal.ii","sex"))
    if (both.sexes && sex == "F") return(c("autosomal.iG","autosomal.iR","autosomal.ii","chrx","chry"))
    if (!both.sexes && sex == "M") return(c("genomic.iG", "genomic.iR", "genomic.ii"))
    if (!both.sexes && sex == "F") return(c("not.y.iG", "not.y.iR", "not.y.ii","chry"))
    stop("invalid input", "sex =", sex, "both.sexes =", both.sexes)
}



extract.controls <- function(rg, probes=meffil.probe.info(), verbose=F) {
    stopifnot(is.rg(rg))

    x.mean <- function(x, na.rm=T) {
        stopifnot(length(x) > 1)
        mean(x,na.rm=na.rm)
    }
    x.which <- function(x) {
        i <- which(x)
        stopifnot(length(i) > 0)
        i
    }

    msg(verbose=verbose)
    probes.G <- probes[x.which(probes$dye == "G"),]
    probes.R <- probes[x.which(probes$dye == "R"),]
    probes.G <- probes.G[match(rownames(rg$G), probes.G$address),]
    probes.R <- probes.R[match(rownames(rg$R), probes.R$address),]

    bisulfite2 <- x.mean(rg$R[x.which(probes.R$target == "BISULFITE CONVERSION II"), "Mean"])

    bisulfite1.G <- rg$G[x.which(probes.G$target == "BISULFITE CONVERSION I"
                               & probes.G$ext
                               %in% sprintf("BS Conversion I%sC%s", c(" ", "-", "-"), 1:3)),"Mean"]
    bisulfite1.R <- rg$R[x.which(probes.R$target == "BISULFITE CONVERSION I"
                               & probes.R$ext %in% sprintf("BS Conversion I-C%s", 4:6)),"Mean"]
    bisulfite1 <- x.mean(bisulfite1.G + bisulfite1.R)

    stain.G <- rg$G[x.which(probes.G$target == "STAINING" & probes.G$ext == "Biotin (High)"),"Mean"]

    stain.R <- rg$R[x.which(probes.R$target == "STAINING" & probes.R$ext == "DNP (High)"),"Mean"]

    extension.R <- rg$R[x.which(probes.R$target == "EXTENSION"
                              & probes.R$ext %in% sprintf("Extension (%s)", c("A", "T"))),"Mean"]
    extension.G <- rg$G[x.which(probes.G$target == "EXTENSION"
                              & probes.G$ext %in% sprintf("Extension (%s)", c("C", "G"))),"Mean"]

    hybe <- rg$G[x.which(probes.G$target == "HYBRIDIZATION"),"Mean"]

    targetrem <- rg$G[x.which(probes.G$target %in% "TARGET REMOVAL"),"Mean"]

    nonpoly.R <- rg$R[x.which(probes.R$target == "NON-POLYMORPHIC"
                            & probes.R$ext %in% sprintf("NP (%s)", c("A", "T"))),"Mean"]

    nonpoly.G <- rg$G[x.which(probes.G$target == "NON-POLYMORPHIC"
                            & probes.G$ext %in% sprintf("NP (%s)", c("C", "G"))),"Mean"]

    spec2.G <- rg$G[x.which(probes.G$target == "SPECIFICITY II"),"Mean"]
    spec2.R <- rg$R[x.which(probes.R$target == "SPECIFICITY II"),"Mean"]
    spec2.ratio <- x.mean(spec2.G,na.rm=T)/x.mean(spec2.R,na.rm=T)

    ext <- sprintf("GT Mismatch %s (PM)", 1:3)
    spec1.G <- rg$G[x.which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext),"Mean"]
    spec1.Rp <- rg$R[x.which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext),"Mean"]
    spec1.ratio1 <- x.mean(spec1.Rp,na.rm=T)/x.mean(spec1.G,na.rm=T)

    ext <- sprintf("GT Mismatch %s (PM)", 4:6)
    spec1.Gp <- rg$G[x.which(probes.G$target == "SPECIFICITY I" & probes.G$ext %in% ext),"Mean"]
    spec1.R <- rg$R[x.which(probes.R$target == "SPECIFICITY I" & probes.R$ext %in% ext),"Mean"]
    spec1.ratio2 <- x.mean(spec1.Gp,na.rm=T)/x.mean(spec1.R,na.rm=T)

    spec1.ratio <- (spec1.ratio1 + spec1.ratio2)/2

    normA <- x.mean(rg$R[x.which(probes.R$target == "NORM_A"),"Mean"], na.rm = TRUE)
    normT <- x.mean(rg$R[x.which(probes.R$target == "NORM_T"),"Mean"], na.rm = TRUE)
    normC <- x.mean(rg$G[x.which(probes.G$target == "NORM_C"),"Mean"], na.rm = TRUE)
    normG <- x.mean(rg$G[x.which(probes.G$target == "NORM_G"),"Mean"], na.rm = TRUE)

    dye.bias <- (normC + normG)/(normA + normT)

    probs <- c(0.01, 0.5, 0.99)
    oob.G <- quantile(rg$G[with(probes.G, x.which(target == "OOB")),"Mean"], na.rm=T, probs=probs)
    oob.R <- quantile(rg$R[with(probes.R, x.which(target == "OOB")),"Mean"], na.rm=T, probs=probs)
    oob.ratio <- oob.G[["50%"]]/oob.R[["50%"]]

    c(bisulfite1=bisulfite1,
      bisulfite2=bisulfite2,
      extension.G=extension.G,
      extension.R=extension.R,
      hybe=hybe,
      stain.G=stain.G,
      stain.R=stain.R,
      nonpoly.G=nonpoly.G,
      nonpoly.R=nonpoly.R,
      targetrem=targetrem,
      spec1.G=spec1.G,
      spec1.R=spec1.R,
      spec2.G=spec2.G,
      spec2.R=spec2.R,
      spec1.ratio1=spec1.ratio1,
      spec1.ratio=spec1.ratio,
      spec2.ratio=spec2.ratio,
      spec1.ratio2=spec1.ratio2,
      normA=normA,
      normC=normC,
      normT=normT,
      normG=normG,
      dye.bias=dye.bias,
      oob.G=oob.G,
      oob.ratio=oob.ratio)
}

identify.bad.probes.beadnum <- function(rg, threshold=3, probes=meffil.probe.info(), verbose=F) {
    msg(verbose=verbose)

    cpg.sites <- unique(na.omit(probes$name))
 
    targets <- c(M="M", U="U")
    n.beads <- sapply(targets, function(target) {
        probes <- probes[which(probes$target == target),]
        probes <- probes[match(cpg.sites, probes$name),]
        r.idx <- which(probes$dye == "R")
        g.idx <- which(probes$dye == "G")        
        n.beads <- rep(NA, nrow(probes))
        n.beads[r.idx] <- rg$R[match(probes$address[r.idx], rownames(rg$R)), "NBeads"]
        n.beads[g.idx] <- rg$G[match(probes$address[g.idx], rownames(rg$G)), "NBeads"]
        n.beads
    })
    n.beads <- pmin(n.beads[,"U"], n.beads[,"M"])
    names(n.beads) <- cpg.sites
    n.beads[which(n.beads < threshold)]
}

identify.bad.probes.detectionp <- function(rg, threshold=0.01, probes=meffil.probe.info(), verbose=F) {
    msg(verbose=verbose)

    dyes <- c(R="R", G="G")
    negative.intensities <- lapply(dyes, function(color) {
        addresses <- with(probes, address[which(dye == color & target == "NEGATIVE")])
        rg[[color]][which(rownames(rg[[color]]) %in% addresses),"Mean"]
    })
    negative.med <- lapply(negative.intensities, median, na.rm=T)
    negative.sd  <- lapply(negative.intensities, mad, na.rm=T)

    cpg.sites <- unique(na.omit(probes$name))

    targets <- c(M="M", U="U")
    detection.stats <- lapply(targets, function(target) {
        probes <- probes[which(probes$target == target),]
        probes <- probes[match(cpg.sites, probes$name),]
        r.idx <- which(probes$dye == "R")
        g.idx <- which(probes$dye == "G")

        intensity <- rep(NA, nrow(probes))
        intensity[r.idx] <- rg$R[match(probes$address[r.idx], rownames(rg$R)),"Mean"]
        intensity[g.idx] <- rg$G[match(probes$address[g.idx], rownames(rg$G)),"Mean"]

        med <- rep(NA, nrow(probes))
        med[r.idx] <- negative.med$R
        med[g.idx] <- negative.med$G

        sd <- rep(NA, nrow(probes))
        sd[r.idx] <- negative.sd$R
        sd[g.idx] <- negative.sd$G

        data.frame(intensity=intensity, med=med, sd=sd)
    })

    p.value <- with(detection.stats, 1-pnorm(M$intensity + U$intensity,
                                             mean=M$med + U$med,
                                             sd=M$sd + U$sd))
    names(p.value) <- cpg.sites
    p.value[which(p.value > threshold)]
}


get.snp.probes <- function(rg, probes=meffil.probe.info(), verbose=F) {
    msg(verbose=verbose)
    probes <- probes[which(probes$target %in% c("UG","MG")),]
    probes$target <- substring(probes$target, 1, 1)
    snp.mu <- rg.to.mu(rg, probes)
    snp.mu$M/(snp.mu$M + snp.mu$U + 100)
}


