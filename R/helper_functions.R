
########## Load environment ##########
env_check <- function(){
    # default packages
    cran_pkg <- data.frame(
        pkg = c("getopt", "bnlearn", "parallel", "dplyr", "yardstick"),
        ver = c("1.20.4", "5.0", "4.4.1", "1.1.4", "1.3.2")
    )
    bioc_pkg <- data.frame(
        pkg = c("GenomicRanges", "rtracklayer", "GenomicAlignments"),
        ver = c("1.56.1", "1.64.0", "1.40.0")
    )
    # check environment
    installed <- rownames(installed.packages())
    ### install CRAN packages
    if (!all(cran_pkg$pkg %in% installed)) {
        if (!requireNamespace("remotes", quietly = T)) {
            suppressMessages(install.packages("remotes"))
        }
        if (!requireNamespace("remotes", quietly = T)) {
            stop("Unable to install package remotes.")
        }
        for(pkg in cran_pkg$pkg) {
            if (!requireNamespace(pkg, quietly = T)) {
                id <- which(pkg %in% cran_pkg$pkg)
                suppressMessages(remotes::install_version(pkg, version = cran_pkg$ver[id]))
                if(!requireNamespace(pkg, quietly = T)){
                    suppressMessages(remotes::install_version(pkg))
                }
            }
        }
    }
    installed <- rownames(installed.packages())
    if(!all(cran_pkg$pkg %in% installed)) {
        stop("Unable to install all required CRAN packages.")
    }
    ### install Bioc packages
    if (!all(bioc_pkg$pkg %in% installed)) {
        if (!requireNamespace("BiocManager", quietly = T)) {
            suppressMessages(install.packages("BiocManager"))
        }
        if (!requireNamespace("BiocManager", quietly = T)) {
            stop("Unable to install package BiocManager.")
        }
        for (pkg in bioc_pkg$pkg) {
            if (!requireNamespace(pkg, quietly = T)) {
                id <- which(pkg %in% bioc_pkg$pkg)
                suppressMessages(remotes::install_version(pkg, version = bioc_pkg$ver[id]))
                if (!requireNamespace(pkg, quietly = T)) {
                    suppressMessages(remotes::install_version(pkg))
                }
            }
        }
    }
    installed <- rownames(installed.packages())
    if(!all(bioc_pkg$pkg %in% installed)){
        stop("Unable to install all required BiocManager packages.")
    }
    # load all required packages
    for(pkg in c(cran_pkg$pkg, bioc_pkg$pkg)) {
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
    return(0)
}

# env_check()


########## Parameter settings ##########
mat <- matrix(c(
    # select functions
    "peakcall", "C", 0, "logical", "run Peak Calling/ Peak Refinement for ATAC-seq",
    "train", "R", 0, "logical", "train RegBay prediction model",
    "predict", "P", 0, "logical", "predict regulatory elements using RegBay pre-trained model",
    # input parameters
    "input", "i", 1, "character", "path to input files (.txt format, separated by comma)",
    "pval", "p", 1, "double", "p-value threshold for peak calling/refinement (DEFAULT:0.05)", 0.05,
    "cutoff", "l", 1, "integer", "cutoff length for broad peaks (DEFAULT:1000)", 1000,
    "window", "w", 1, "integer", "window size of peak calling (DEFAULT:10000)", 10000,
    "flank", "f", 1, "integer", "extend length to search for histone modification signals (DEFAULT: 500)", 500,
    # pre-trained model
    "model", "m", 1, "character", "path to pre-train RegBay model (.rds format)",
    # general parameters
    "cores", "t", 1, "integer", "number of cores to use (DEFAULT:1)", 1,
    "verbose", "v", 1, "logical", "show detailed message", TRUE,
    "outdir", "o", 1, "character", "output directory (DEFAULT:current folder)", getwd(),
    "help", "h", 0, "logical", "show help message"
), ncol = 5, byrow = T)

# getopt(mat)

##### Check files #####
check_file <- function(fn) {
    # check single file
    if (!file.exists(fn)) {
        stop("File does not exist.")
    }
    if (file.info(fn)$size == 0) {
        stop("Empty file is detected.")
    }
    return(0)
}

check_dir <- function(dir) {
    if (!dir.exists(dir)) {
        stop("Directory does not exist.")
    }
    return(0)
}

