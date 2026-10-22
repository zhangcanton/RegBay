cat(
    "\n
    ============================================================
                         RegBay v1.0 (2026)
       Regulatory element prediction based on Bayesian network
    ============================================================\n\n
    "
)

########## path to project resources ##########
get_source_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    # Check whether script path exists
    if (length(file_arg) == 0) {
        stop(
            "Cannot determine script path. This function must be run using Rscript or R --file."
        )
    } else if (length(file_arg) > 1) {
        # Handle unexpected multiple matches
        stop("Multiple --file arguments detected.")
    }
    # Extract full path to script
    script_path <- tryCatch(
        normalizePath(sub("^--file=", "", file_arg)),
        error = function(e) {
            stop("Script file does not exist or path is invalid.")
        }
    )
    # Path to project directory
    proj_dir <- dirname(script_path)
    return(proj_dir)
}

path <- get_source_dir()
source(paste0(path, "/R/helper_functions.R"))

########## user-input parameters ##########
if (!requireNamespace("getopt", quietly = T)) {
    suppressMessages(install.packages("getopt"))
}
suppressPackageStartupMessages(library(getopt))

opt <- getopt(mat)
if (all(sapply(c(opt$help, opt$peakcall, opt$train, opt$predict), is.null))) {
    stop("Must specify one function module (peakCall, train, predict) or \'help\'")
}

### Help messages
if (!is.null(opt$help)) {
    cat(
        "
    Usage: Rscript regBay.R [options]

    Options:

      Module:
        -C, --peakcall  Peak calling/refining module
        -R, --train     Training module for RegBay regulatory element model
        -P, --predict   Prediction module for regulatory element using RegBay

      I/O option:
        -i, --input     Path to input files (.txt format)
        -m, --model     Path to RegBay model (.rds format)
        -o, --outdir    Output dircetory (DEFAULT:current folder)

      Algorithm parameter:
        -p, --pval      P-value for peak calling (DEFAULT:0.05)
        -l, --cutoff    Cutoff length for broad peaks (DEFAULT:1000)
        -w, --window    Window size of peak calling (DEFAULT:10000)
        -f, --flank     Size of flanking regions (DEFAULT:500)

      Performance option:
        -t, --cores     Number of cores to use (DEFAULT:1)

      Miscellaneous:
        -h, --help      Show help message
        -v, --verbose   Show detailed message\n
"
    )
    q(status = 0)
}

########## set up required environment ##########
flag <- env_check()
source(paste0(path, "/R/functions.R"))

########## initiate funciton modules ##########
# Peak calling/Peak refinement
if (!is.null(opt$peakcall)) {
    opt <- set_default(opt)
    source(paste0(path, "/R/peak_calling.R"))
}

# Train RegBay model
if (!is.null(opt$train)) {
    opt <- set_default(opt)
    source(paste0(path, "/R/regbay_train.R"))
}

# Predict regulatory elements using RegBay model
if (!is.null(opt$predict)) {
    opt <- set_default(opt)
    source(paste0(path, "/R/regbay_predict.R"))
}
