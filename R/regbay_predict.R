# ============================================================
# Predict regulatory elements using RegBay
# ============================================================

cat(
    "
    +++++  Predict regulatory elements using pre-trained model   +++++
    "
)

### Start time
s.time <- Sys.time()

### check output dir
flag <- check_dir(opt$outdir)

### Read input meta txt
meta <- read_meta(opt)

### Sequencing BAM files (mandatory)
if (opt$verbose) {
    cat(
    "
    >>> Counting sequencing coverage.")
}
# ATAC-seq coverage
raws <- seq_cov(meta$atac)
raws <- raws[sapply(raws, function(x) length(x@values)) > 1] # remove empty chr
# HM coverage
hms <- seq_cov(meta, cores = opt$cores)
if(!all(tolower(c("H3K27ac", "H3K4me1", "H3K4me3")) %in% names(hms))) {
    stop("Missing histone modification profiles.")
}
if(!"input" %in% names(hms)) {
    cat(
    "
    [INFO] No input/control for histone modifications is detected."
    )
    if(opt$verbose){
        cat(
    "
    >>> Running analysis without input/control.")
    }
} else{
    cat(
    "
    [INFO] Found input/control for histone modifications."
    )
    if (opt$verbose) {
        cat(
    "
    >>> Running analysis with input/control."
        )
    }
}

### Normalization
if (opt$verbose) {
    cat(
    "
    >>> Normalizing sequencing counts.")
}
# library size
tot <- mclapply(meta[names(hms)], function(bam){
  cots <- idxstatsBam(bam)
  return(sum(cots$mapped))
}, mc.cores = opt$cores)
tot <- unlist(tot)
# Log2 normalization
if (any(names(hms) %in% "input")) {
  # with input/none
  names(tot) <- names(hms)
  scales <- unlist(tot[which(names(tot) != 'input')]) / tot["input"]
  norms <- mclapply(seq_along(scales), function(i){
    # log2 normalization to input signals
    inputs <- hms[['input']] * scales[i]
    counts <- log2((hms[[i]] + 1) / (hms[["input"]] + 1))
    # range normalization (skip potential outliers)
    upper <- round(max(quantile(counts, 0.995))) + 1
    lower <- round(min(quantile(counts, 0.005))) - 1
    counts[counts > 0] <- counts[counts > 0] / upper * 20
    counts[counts > 20] <- 20
    counts[counts < 0] <- counts[counts < 0] / lower * -10
    counts[counts < -10] <- -10
    counts <- counts[sapply(counts, function(x) length(x@values)) > 1] # in case some chr are empty, kick out those
    return(counts)
  }, mc.cores = opt$cores)
  names(norms) <- names(scales)
  cov <- log2(raws + 1)
  scales <- round(max(quantile(cov, 0.995))) + 1
  cov <- cov / scales * 20
  cov[cov > 20] <- 20
  norms <- append(norms, list(atac = cov))
} else {
  # without input/none
  norms <- mclapply(names(hms), function(i) {
    counts <- log2(hms[[i]] + 1)
    upper <- round(max(quantile(counts, 0.995))) + 1
    counts <- counts / upper * 20
    counts[counts > 20] <- 20
    counts <- counts[sapply(counts, function(x) length(x@values)) > 1] # in case some chr are empty, kick out those
    return(counts)
  }, mc.cores = cores)
  names(norm_counts) <- names(hms)
  cov <- log2(raws + 1)
  scales <- round(max(quantile(counts, 0.995))) + 1
  cov <- cov / scales * 20
  cov[cov > 20] <- 20
  norms <- append(norms, list(atac = cov))
}

### call ATAC-seq peaks
peaks <- tryCatch(
    import(paste0(opt$outdir, "/peaks.bed")),
    error = function(e) {
        stop("Unable to read peak file.")
    }
)

### flanking regions
prev <- resize(peaks, width = width(peaks) + opt$flank, fix = "end")
prev <- resize(prev, width = opt$flank, fix = "start")
folw <- resize(peaks, width = width(peaks) + opt$flank, fix = "start")
folw <- resize(folw, width = opt$flank, fix = "end")
overs <- IRanges::intersect(prev, peaks)
prev <- IRanges::setdiff(prev, overs)
overs <- IRanges::intersect(folw, peaks)
folw <- IRanges::setdiff(folw, overs)
### coverage at peaks
prange <- split(ranges(peaks), peaks@seqnames)
pnorms <- lapply(names(norms), function(n){
  mclapply(names(prange), function(i){
    deps <- extractList(norms[[n]][[i]], prange[[i]])
    deps <- max(deps)
    return(deps)
  }, mc.cores = opt$cores)
})
pnorms <- lapply(pnorms, unlist)
pnorms <- data.frame(matrix(unlist(pnorms), ncol = 4, byrow = F))
names(pnorms) <- names(norms)
pprev <- split(ranges(prev), peaks@seqnames)
nprev <- lapply(names(norms)[names(norms) != 'atac'], function(n) {
  mclapply(names(pprev), function(i) {
    deps <- extractList(norms[[n]][[i]], pprev[[i]])
    deps <- max(deps)
    return(deps)
  }, mc.cores = opt$cores)
})
nprev <- lapply(nprev, unlist)
nprev <- data.frame(matrix(unlist(nprev), ncol = 3, byrow = F))
names(nprev) <- names(norms)[names(norms) != "atac"]
pfolw <- split(ranges(folw), peaks@seqnames)
nfolw <- lapply(names(norms)[names(norms) != 'atac'], function(n) {
  mclapply(names(pfolw), function(i) {
    deps <- extractList(norms[[n]][[i]], pfolw[[i]])
    deps <- max(deps)
    return(deps)
  }, mc.cores = opt$cores)
})
nfolw <- lapply(nfolw, unlist)
nfolw <- data.frame(matrix(unlist(nfolw), ncol = 3, byrow = F))
names(nfolw) <- names(norms)[names(norms) != "atac"]
neigh <- pmax(nprev, nfolw)
pnorms[, 1:3] <- pmax(neigh, pnorms[, 1:3])
pnorms$me_ratio <- pnorms$h3k4me3 - pnorms$h3k4me1
pnorms$me_ratio <- cut(pnorms$me_ratio, breaks = c(-Inf, -1, 1, Inf), labels = c(0, 1, 2))

### map sequencing profiles to candidate genomic regions
candidates <- peaks
mcols(candidates) <- pnorms
# factorize counts
mcols(candidates)[, -grep("actP|actE|me_ratio", names(mcols(candidates)))] <- lapply(
  mcols(candidates)[, -grep("actP|actE|me_ratio", names(mcols(candidates)))],
  function(x) as.integer(round(x))
)
candidates$atac <- factor(candidates$atac, levels = seq(0, 20))
mcols(candidates)[, -grep("actP|actE|me_ratio|ATAC", names(mcols(candidates)))] <- lapply(
  mcols(candidates)[, -grep("actP|actE|me_ratio|ATAC", names(mcols(candidates)))],
  function(x) {
    tmp <- order(x)
    x <- factor(sort(x), levels = seq(-10, 20))
    x <- x[order(tmp)]
  }
)
# transform to data.frame
candis <- data.frame(mcols(candidates))


### Bayesian Network
if (opt$verbose) {
    cat(
    "
    >>> Predicting regulatory elements using Bayesian network model."
    )
}
# pre-trained model
fit.bn <- tryCatch(
    readRDS(opt$model),
    error = function(e) {
        stop("Unable to load prediction model.")
    }
)
# predict regulatory elements on candidates
pActE <- predict(fit.bn, node = "actE", data = candis, prob = T, method = "parents")
pActP <- predict(fit.bn, node = "actP", data = candis, prob = T, method = "parents")

### output genomic regions with predicted probabilities
# map and output to BigWig files
gr <- candidates
mcols(gr) <- NULL
gr$score <- attr(pActP, "prob")[2, ]
export.bedGraph(gr[gr$score > 0.5], paste0(opt$outdir, "/predict.promoter.bedGraph"))
gr$score <- attr(pActE, "prob")[2, ]
export.bedGraph(gr[gr$score > 0.5], paste0(opt$outdir, "/predict.enhancer.bedGraph"))
# output to RDS file
mcols(gr) <- NULL
gr$pmo <- attr(pActP, "prob")[2, ]
gr$enc <- attr(pActE, "prob")[2, ]
saveRDS(gr, paste0(opt$outdir, "/predictions.RDS"))
# closing message
cat(
    "
    DONE: Predicting regulatory elements. Run time:", round(difftime(Sys.time(), s.time, units = "mins"), digits = 2), "mins\n"
)
