#' Regulatory element prediction based on Bayesian network
#' @description
#' RegBay predicts regulatory element using ATAC-seq and histone modification ChIP-seq profiles.
#' @param fn_list variable names in 1st column and path to files in 2nd column, separated by comma (no space)
#' files include path to bam files, annotation bed files, and the optional ATAC peak file
#' @param genome genome name or genome files
#' @param p_val p-value threshold
#' @param cores number of available cores, for parallel running
#' @return prediction model
#' @examples
#' bn_model <- bn_train(meta = data.meta, genome = "hg38", p_val = 0.05, cores = 8)
#' @export
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import BSgenome.Mmusculus.UCSC.mm10
#' @import GenomicRanges
#' @import rtracklayer
#' @import bamsignals
#' @import GenomicAlignments
#' @import preprocessCore
#' @importFrom parallel mclappy
#' @import bnlearn
#' @import dplyr


##### Part 1: Data Preparation #####
### parameters
# opt

### Sequencing BAM files (mandatory)
# file check
bam_cond <- c("ATAC", "H3K27ac", "H3K4me1", "H3K4me3", "input")
bam_cond <- tolower(bam_cond) # NOT tested
# peak_fn <- "/project/EnhancerAnnotation/new_files/manuscript/tracks/k562_atac_p05_genrich.narrowPeak"
if (length(bam_cond) < 4) {
  stop("Missing BAM input")
}

### Annotations (mandatory)
# file check




##### Part 2: ATAC-seq processing and peak calling #####
print('Start processing ...')
print('>>> Calling ATAC-seq peaks ...')
##### ATAC-seq coverage
raws <- readGAlignments(bam_fn[which(bam_cond == 'atac')])
cov <- coverage(raws)
cov <- cov[grep("^chr[0-9,X]{1,2}$", names(cov))]
cov <- cov[sapply(cov, function(x) length(x@values)) > 1] # in case some chr are empty, kick out those
### call ATAC-seq peaks
path <- dirname(getwd())
source(paste0(path,'/R/peak_calling.R'))
if(){
  peaks <- peak_calling(cov, peak_fn = peak_fn, p_val = p_val, cores = cores)
} else {
  peaks <- peak_calling(cov, p_val = p_val, cores = cores)
}


                  
##### Part 3: Coverage at BAM files #####
print(">>> Counting read depth from BAM files ...")
### raw coverage from BAM files
raws <- cov
hms <- mclapply(bam_fn[which(bam_cond != "atac")], function(bam){
  cov <- coverage(bam)
  cov <- cov[grep("^chr[0-9,X]{1,2}$", names(cov))]
  return(cov)
}, mc.cores = cores)
names(hms) <- bam_cond[which(bam_cond != 'atac')]

tot <- mclapply(bam_fn[which(bam_cond != "atac")], function(bam){
  cots <- idxstatsBam(bam)
  return(sum(cots$mapped))
}, mc.cores = cores)
tot <- unlist(tot)

### read depth normalization
# Log2 normalization
if (any(bam_cond %in% "input")) {
  # with input/none
  names(tot) <- bam_cond[which(bam_cond != 'atac')]
  scales <- unlist(tot[which(names(tot) != 'input')]) / tot[["input"]]
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
  }, mc.cores = cores)
  names(norms) <- names(scales)
  cov <- log2(raws + 1)
  scales <- round(max(quantile(cov, 0.995))) + 1
  cov <- cov / scales * 20
  cov[cov > 20] <- 20
  norms <- append(norms, list(atac = cov))
} else {
  # without input/none
  norms <- mclapply(bam_cond, function(i) {
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



##### Part 4: Construction of feature matrix #####
### flanking regions
# Note: extend 500 bp or reach another ATAC-seq peak
prev <- resize(peaks, width = width(peaks) + 500, fix = "end")
prev <- resize(prev, width = 500, fix = "start")
folw <- resize(peaks, width = width(peaks) + 500, fix = "start")
folw <- resize(folw, width = 500, fix = "end")
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
  }, mc.cores = cores)
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
  }, mc.cores = cores)
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
  }, mc.cores = cores)
})
nfolw <- lapply(nfolw, unlist)
nfolw <- data.frame(matrix(unlist(nfolw), ncol = 3, byrow = F))
names(nfolw) <- names(norms)[names(norms) != "atac"]
neigh <- pmax(nprev, nfolw)
pnorms[, 1:3] <- pmax(neigh, pnorms[, 1:3])
pnorms$me_ratio <- pnorms$h3k4me3 - pnorms$h3k4me1
pnorms$me_ratio <- cut(pnorms$me_ratio, breaks = c(-Inf, -1, 1, Inf), labels = c(0, 1, 2))

### labelling peaks
candidates <- peaks
mcols(candidates) <- pnorms
actE_index <- candidates %over% active_enc
inaE_index <- candidates %over% inactive_enc & candidates %outside% active_enc
actP_index <- candidates %over% active_pmo
inaP_index <- candidates %over% inactive_pmo & candidates %outside% active_pmo

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
mcols(candidates)[, -grep("actP|actE|me_ratio|ATAC", names(mcols(candidates)))] <- lapply(
  mcols(candidates)[, -grep("actP|actE|me_ratio|ATAC", names(mcols(candidates)))],
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
### transform to data.frame
cand.df <- data.frame(mcols(candidates))
cand.df <- cand.df[which(!is.na(cand.df$actE) & !is.na(cand.df$actP)), ]

### Splitting training and test data
set.seed(24)
actE <- which(cand.df$actE == 1)
actP <- which(cand.df$actP == 1)
inactive <- which(!(1:nrow(cand.df) %in% unique(actE, actP)))
training <- actE[sample(1:length(actE), round(0.7*length(actE)))]
training <- c(training, actP[sample(1:length(actP), round(0.7*length(actP)))])
training <- c(training, inactive[sample(1:length(inactive), round(0.7*length(inactive)))])
training_ind <- unique(training)
test_ind <- !(1:nrow(cand.df) %in% training_ind)
training <- cand.df[training_ind,]
test <- cand.df[test_ind,]

### Bayesian Network I: structure learning
# Structure 1: with ATAC and active nodes
BN_graph <- empty.graph(nodes = c(
  "atac", "h3k27ac", "h3k4me1", "h3k4me3",
  "actP", "actE", "me_ratio"
))
arc.set <- matrix(
  c(
    "h3k4me1", "actP", "h3k4me3", "actP",
    "h3k4me1", "actE", "h3k4me3", "actE",
    "atac", "actP", "atac", "actE",
    "h3k4me1", "me_ratio", "h3k4me3", "me_ratio",
    "h3k27ac", 'actP', 'h3k27ac', 'actE',
    "me_ratio", "actP", "me_ratio", "actE"
  ),
  byrow = TRUE, ncol = 2,
  dimnames = list(NULL, c("from", "to"))
)
arcs(BN_graph) <- arc.set

### Bayesian Network II: parameter learning
fit.bn <- bn.fit(BN_graph, data = training, method = "bayes")

### Bayesian Network III: prediction and validation
## Predict nodes
pActE <- predict(fit.bn, node = "actE", data = test, prob = T, method = "parents")
pActP <- predict(fit.bn, node = "actP", data = test, prob = T, method = "parents")

### Accuracy
acc.label <- list(enc = NULL, pmo = NULL)
tab <- table(test$actE, pActE)
acc.label[["pmo"]] <- paste("Pmo (Acc = ", toString(round((sum(diag(tab)) / sum(tab)) * 100, 2)), ")", sep = "")
tab <- table(test$actP, pActP)
acc.label[["enc"]] <- paste("Enc (Acc = ", toString(round((sum(diag(tab)) / sum(tab)) * 100, 2)), ")", sep = "")

### PR curve
pr.label <- list(enc = NULL, pmo = NULL)
pred <- data.frame(pAct = pActP, tAct = test$actP, prob = attr(pActP, "prob")[2, ])
pred$pAct <- factor(pred$pAct, level = c(1, 0))
pred$tAct <- factor(pred$tAct, level = c(1, 0))
pred2 <- data.frame(pAct = pActE, tAct = test$actE, prob = attr(pActE, "prob")[2, ])
pred2$pAct <- factor(pred2$pAct, level = c(1, 0))
pred2$tAct <- factor(pred2$tAct, level = c(1, 0))
auc_pred <- pr_auc(data = pred, truth = tAct, prob)$.estimate
pr.label[["pmo"]] <- paste("Pmo (AUPRC = ", toString(round(auc_pred, 4)), ")", sep = "")
auc2_pred <- pr_auc(data = pred2, truth = tAct, prob)$.estimate
pr.label[["enc"]] <- paste("Enc (AUPRC = ", toString(round(auc2_pred, 4)), ")", sep = "")

### ROC curve
roc.label <- list(enc = NULL, pmo = NULL)
# pred <- data.frame(pAct = pActP, tAct = test$actP, prob = attr(pActP, "prob")[2, ])
# pred2 <- data.frame(pAct = pActE, tAct = test$actE, prob = attr(pActE, "prob")[2, ])
auc_pred <- roc_auc(data = pred, tAct, prob)
auc_pred <- auc_pred$.estimate
roc.label[["pmo"]] <- paste("Pmo (AUC=", toString(round(1 - auc_pred, 4)), ")", sep = "")
auc2_pred <- roc_auc(data = pred2, tAct, prob)
auc2_pred <- auc2_pred$.estimate
roc.label[["enc"]] <- paste("Enc (AUC=", toString(round(1 - auc2_pred, 4)), ")", sep = "")

# output
print("Accuracy:")
do.call(cat, c(acc.label, sep = "\n"))
print("AUPRC:")
do.call(cat, c(pr.label, sep = "\n"))
print("AUC:")
do.call(cat, c(roc.label, sep = "\n")

### output prediction model
return(fti.bn)
