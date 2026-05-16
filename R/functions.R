
##### Read input file #####
read_meta <- function(opt){
    # check input file
    if (is.null(opt$input)) {
        stop("Augment \"input\" is missing.")
    }
    meta <- tryCatch(
        read.table(opt$input, sep = ","),
        error = function(e) {
            stop("Unable to read input txt file.")
        }
    )
    colnames(meta) <- c("var", "fn")
    meta$var <- tolower(meta$var)
    
    # check ATAC
    if (!("atac" %in% meta$var)) {
        stop("Missing BAM file for ATAC-seq.")
    }
    flag <- check_file(meta$fn[which("atac" %in% meta$var)])
    # if ("peak_fn" %in% meta$var) {
    #     flag <- check_file(meta$fn[which("peak_fn" %in% meta$var)])
    #     print("Detect a ATAC-seq peak file. Skip peak calling.")
    # }
    
    # check HMs
    if(!is.null(opt$train) | !is.null(opt$predict)){
        if (!(("h3k27ac" %in% meta$var) & ("h3k4me3" %in% meta$var) & ("h3k4me1" %in% meta$var))) {
            stop("Missing BAM file for ChIP-seq.")
        }
        if (!("input" %in% meta$var)) {
            stop("Missing control file for ChIP-seq.")
        }
        flag <- lapply(meta$fn[which(meta$var %in% c('h3k27ac', 'h3k4me1', 'h3k4me3', 'input'))], check_file)
    }

    # check annotations
    if (!is.null(opt$train) | !is.null(opt$predict)) {
        if (!(("active_pmo" %in% meta$var) & ("inactive_pmo" %in% meta$var) & ("active_enc" %in% meta$var) & ("inactive_enc" %in% meta$var))) {
            stop("Missing BED file for annotations.")
        }
        flag <- lapply(meta$fn[which(meta$var %in% c("active_pmo", "inactive_pmo", "active_enc", "inactive_enc"))], check_file)
    }

    # set variables
    fns <- setNames(as.list(meta$fn), meta$var)
    return(fns)
}


##### Build annotation sets #####
get_anno <- function(flist){
    flist <- flist[grepl("act", names(flist))]
    if(length(flist) < 4 | !all(
        c('active_pmo', 'inactive_pmo', 'active_enc', 'inactive_enc') %in% names(flist))) {
        stop("Missing annotation sets.")
    }
    gr <- mapply(
        function(v, f) assign(v, import.bed(f)),
        names(flist),
        flist
    )
    return(gr)
}


##### Count sequence coverage #####
seq_cov <- function(fn, cores = 1) {
    if (is.character(fn)) {
        cov <- coverage(fn)
        cov <- cov[grep("^chr[0-9,X]{1,2}$", names(cov))]
    } else if (is.list(fn)) {
        ids <- sapply(c("h3k27ac", "h3k4me1", "h3k4me3", "input"), function(ch) grep(ch, names(fn)))
        cov <- mclapply(
            fn[ids[sapply(ids, length) != 0]],
            function(f) {
                dep <- coverage(f)
                dep <- dep[grep("^chr[0-9,X]{1,2}$", names(dep))]
                return(dep)
            },
            mc.cores = cores
        )
        names(cov) <- c("h3k27ac", "h3k4me1", "h3k4me3", "input")[sapply(ids, length) != 0]
    }
    return(cov)
}


### Bayesian network structure
bn_net <- function() {
    netw <- empty.graph(nodes = c(
        "atac", "h3k27ac", "h3k4me1", "h3k4me3",
        "actP", "actE", "me_ratio"
    ))
    arc.set <- matrix(
        c(
            "h3k4me1", "actP", "h3k4me3", "actP",
            "h3k4me1", "actE", "h3k4me3", "actE",
            "atac", "actP", "atac", "actE",
            "h3k4me1", "me_ratio", "h3k4me3", "me_ratio",
            "h3k27ac", "actP", "h3k27ac", "actE",
            "me_ratio", "actP", "me_ratio", "actE"
        ),
        byrow = TRUE, ncol = 2,
        dimnames = list(NULL, c("from", "to"))
    )
    arcs(netw) <- arc.set
    return(netw)
}
