#' Regulatory element prediction based on Bayesian network
#' @description
#' RegBay predicts regulatory element using ATAC-seq and histone modification ChIP-seq
#' profiles.
#' @param bam_fn path to bam files including ATAC, H3K27ac, H3K4me1, and H3K4me3
#' @param active_pmo path to a bed file
#' @param inactive_pmo path to a bed file
#' @param active_enc path to a bed file
#' @param inactive_enc path to a bed file
#' @param peak_fn optional, address to peak files, BED or narrowPeak
#' @param p_val p-value threshold
#' @param cores number of available cores, for parallel running
#' @return BED or BigWig file genomic regions with predicted probabilities of promoters and enhancers
#' @examples
#' bn_train()
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
