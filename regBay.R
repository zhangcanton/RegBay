
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


########## set up required environment ##########
env_check()

########## user-input parameters ##########
opt <- getopt(mat)
if(all(sapply(c(opt$help, opt$peakCall, opt$train, opt$predict), is.null))){
  stop("Must specify one function module (peakCall, train, predict) or \'help\'")
}

### Help messages


########## start up funciton modules ##########
# Peak calling/Peak refinement
if (!is.null(opt$peakCall)) source(paste0(path, "/BIN/peak_calling.R"))

# Train RegBay model
if (!is.null(opt$train)) source(paste0(path, "/BIN/regbay_train.R"))

# Predict regulatory elements using RegBay model
if (!is.null(opt$predict)) source(paste0(path, "/BIN/regbay_predict.R"))
