
cat(
    "\n
    ============================================================
                            RegBay v1.0
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
source(paste0(path, "/BIN/basic_functions.R"))
source(paste0(path, "/BIN/functions.R"))


########## set up required environment ##########
env_check()

########## user-input parameters ##########
opt <- getopt(mat)
if(all(sapply(c(opt$help, opt$peakCall, opt$train, opt$predict), is.null))){
  stop("Must specify one function module (peakCall, train, predict) or \'help\'")
}

### Help messages
if(!is.null(opt$help)){
    cat(
        "
        Usage: Rscript regBay.R [options]
        
        Options:
            -h, --help  Show help message
            -o, --outdir    output dircetory (DEFAULT:current folder)
            -t, --cores     number of cores to use (DEFAULT:1)
            -C, --peakCall  Peak calling/refining module
            -R, --train     Train RegBay regulatory element model
            -P, --predict   Predict regulatory element using RegBay
            -i, --input     path to input files (.txt format)
            -p, --pval      p-value for peak calling (DEFAULT:0.05)
            -m, --model     path to RegBay model (.rds format)
            -l, --cutoff    cutoff length for broad peaks (DEFAULT:1000)
            -w, --window    window size of peak calling (DEFAULT:10000)
            -f, --flank     size of flanking regions (DEFAULT:500)
        "
    )
}

########## initiate funciton modules ##########
# Peak calling/Peak refinement
if (!is.null(opt$peakCall)) source(paste0(path, "/BIN/peak_calling.R"))

# Train RegBay model
if (!is.null(opt$train)) source(paste0(path, "/BIN/regbay_train.R"))

# Predict regulatory elements using RegBay model
if (!is.null(opt$predict)) source(paste0(path, "/BIN/regbay_predict.R"))
