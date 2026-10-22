# ============================================================
# Train RegBay prediction model
# ============================================================

cat(
    "
    +++++  Train prediction model for regulatory element   +++++
    "
)

### Start time
s.time <- Sys.time()

### check output dir
flag <- check_dir(opt$outdir)

### Read input meta txt
meta <- read_meta(opt)

### Annotations (mandatory)
if(opt$verbose){
    cat(
    "
    >>> Loading regulatory element annotation sets.")
}
annos <- get_anno(meta)

### Sequencing BAM files (mandatory)
if (opt$verbose) {
    cat(
    "
    >>> Counting sequencing coverage.")
}
# bam_cond <- c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "input")
# bam_cond <- tolower(bam_cond)
# ATAC-seq coverage
raws <- seq_cov(meta$atac)
raws <- raws[sapply(raws, function(x) length(x@values)) > 1] # remove empty chr
# HM coverage
hms <- seq_cov(meta, cores = opt$cores)
if(!all(tolower(c("H3K27ac", "H3K4me1", "H3K4me3")) %in% names(hms))) {
    stop("Missing histone modification profiles.")
}
if (!"input" %in% names(hms)) {
  cat(
    "
    [INFO] No input/control for histone modifications is detected.")
  if (opt$verbose) {
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
    counts <- log2((hms[[i]] + 1) / (inputs + 1))
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
if (is.null(meta$peak_fn)) {
  stop("Missing peak file.")
}
flag <- check_file(meta$peak_fn)
peaks <- tryCatch(
  import(meta$peak_fn),
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

### labelling candidate genomic regions
candidates <- peaks
mcols(candidates) <- pnorms
actE_index <- candidates %over% annos$active_enc
inaE_index <- candidates %over% annos$inactive_enc & candidates %outside% annos$active_enc
actP_index <- candidates %over% annos$active_pmo
inaP_index <- candidates %over% annos$inactive_pmo & candidates %outside% annos$active_pmo

mcols(candidates)["actE"] <- NA
mcols(candidates)["actP"] <- NA
mcols(candidates)[which(actE_index == 1), "actE"] <- 1
mcols(candidates)[which(inaE_index == 1), "actE"] <- 0
mcols(candidates)[which(actP_index == 1), "actP"] <- 1
mcols(candidates)[which(inaP_index == 1), "actP"] <- 0
# in principle, actP and actE should not appear in the same candidate regions
mcols(candidates)[which(mcols(candidates)[, "actP"] == 1 & is.na(mcols(candidates)["actE"])), "actE"] <- 0
mcols(candidates)[which(mcols(candidates)[, "actE"] == 1 & is.na(mcols(candidates)["actP"])), "actP"] <- 0
candidates$actE[which(candidates$actP == 1)] <- 0
# factorize counts
mcols(candidates)[, -grep("actP|actE|me_ratio", names(mcols(candidates)))] <- lapply(
  mcols(candidates)[, -grep("actP|actE|me_ratio", names(mcols(candidates)))],
  function(x) as.integer(round(x))
)
candidates$atac <- factor(candidates$atac, levels = seq(0, 20))
mcols(candidates)[, -grep("actP|actE|me_ratio|atac", names(mcols(candidates)))] <- lapply(
  mcols(candidates)[, -grep("actP|actE|me_ratio|atac", names(mcols(candidates)))],
  function(x) {
    tmp <- order(x)
    x <- factor(sort(x), levels = seq(-10, 20))
    x <- x[order(tmp)]
  }
)
mcols(candidates)[, grep("actP|actE", names(mcols(candidates)))] <- lapply(
  mcols(candidates)[, grep("actP|actE", names(mcols(candidates)))],
  function(x) {
    tmp <- order(x)
    x <- factor(sort(x), levels = c(0, 1))
    x <- x[order(tmp)]
  }
)

### Splitting training and test data
# transform to data.frame
cand.df <- data.frame(mcols(candidates))
cand.df <- cand.df[which(!is.na(cand.df$actE) & !is.na(cand.df$actP)), ]
# divide training and test sets
set.seed(24)
actE <- which(cand.df$actE == 1)
actP <- which(cand.df$actP == 1)
inactive <- which(!(1:nrow(cand.df) %in% unique(c(actE, actP))))
training <- actE[sample(1:length(actE), round(0.7*length(actE)))]
training <- c(training, actP[sample(1:length(actP), round(0.7*length(actP)))])
training <- c(training, inactive[sample(1:length(inactive), round(0.7*length(inactive)))])
training_ind <- unique(training)
test_ind <- !(1:nrow(cand.df) %in% training_ind)
training <- cand.df[training_ind,]
test <- cand.df[test_ind,]

### Bayesian Network
if (opt$verbose) {
    cat(
    "
    >>> Training prediction model based on Bayesian network."
    )
}
# structure
net_graph <- bn_net()
# parameter
fit.bn <- bn.fit(net_graph, data = training, method = "bayes")
# predict nodes
pActE <- predict(fit.bn, node = "actE", data = test, prob = T, method = "parents")
pActP <- predict(fit.bn, node = "actP", data = test, prob = T, method = "parents")

### Metrices
df <- data.frame(enc = rep(0,3), pmo = rep(0,3))
rownames(df) <- c("Accuracy", "AUROC", "AUPRC")
pred <- data.frame(pAct = pActP, tAct = test$actP, prob = attr(pActP, "prob")[2, ])
pred$pAct <- factor(pred$pAct, level = c(1, 0))
pred$tAct <- factor(pred$tAct, level = c(1, 0))
pred2 <- data.frame(pAct = pActE, tAct = test$actE, prob = attr(pActE, "prob")[2, ])
pred2$pAct <- factor(pred2$pAct, level = c(1, 0))
pred2$tAct <- factor(pred2$tAct, level = c(1, 0))
# Accuracy
tab <- table(test$actE, pActE)
df$enc[1] <- round((sum(diag(tab)) / sum(tab)) * 100, 2)
tab <- table(test$actP, pActP)
df$pmo[1] <- round((sum(diag(tab)) / sum(tab)) * 100, 2)
# ROC curve
auc_pred <- roc_auc(data = pred, tAct, prob)$.estimate
df$pmo[2] <- round(auc_pred * 100, 2)
auc2_pred <- roc_auc(data = pred2, tAct, prob)$.estimate
df$enc[2] <- round(auc2_pred * 100, 2)
# PR curve
auc_pred <- pr_auc(data = pred, truth = tAct, prob)$.estimate
df$pmo[3] <- round(auc_pred * 100, 2)
auc2_pred <- pr_auc(data = pred2, truth = tAct, prob)$.estimate
df$enc[3] <- round(auc2_pred * 100, 2)

# output
cat(
        "\n
        *** Performance on test set ***")
cat(sprintf(
  "
          %-26s", "--------------------------"))
cat(sprintf(
        "
          %-8s %-8s %-8s", "Metrics", "Enhancer", "Promoter"))
cat(sprintf(
        "
          %-26s", "--------------------------"))
for (i in rownames(df)) {
  cat(sprintf(
    "
          %-8s %-8s %-8s", i, df[i, "enc"], df[i, "pmo"]
  ))
}
cat(sprintf(
  "
          %-26s", "--------------------------\n"))

### output trained model
saveRDS(fit.bn, paste0(opt$outdir, "/model.RDS"))
cat(
    "
    DONE: Model training. Run time:", round(difftime(Sys.time(), s.time, units = "mins"), digits = 2), "mins\n"
)
