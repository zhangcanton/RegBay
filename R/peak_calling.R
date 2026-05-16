# ============================================================
# ATAC-seq peak calling/ peak refinement
# ============================================================

cat(
    "
    +++++              ATAC-seq peak calling               +++++\n
    "
)

### Start time
s.time <- Sys.time()

### check output dir
flag <- check_dir(opt$outdir)

### Read input meta txt
meta <- read_meta(opt)

### define necessary parameters
if ("peak_fn" %in% names(meta)) {
    flag <- check_file(meta$peak_fn)
    suffix <- sub(".*\\.", "", meta$peak_fn)
    if (suffix != "bed" | suffix != "narrowPeak") {
        stop("Unable to read peak file.")
    }
    if (opt$verbose) {
        cat("Detect a ATAC-seq peak file. Skip peak calling.\n")
    }
    peak_fn <- meta$peak_fn
} else {
    peak_fn <- NULL
    if (opt$verbose) {
        cat("No peak file detected. Start peak calling.\n")
    }
}

### Count ATAC-seq coverage
if (opt$verbose) {
    cat(">>> Count ATAC-seq coverage ...\n")
}
cov <- seq_cov(meta$atac)
cov <- cov[grep("^chr[0-9,X]{1,2}$", names(cov))]
cov <- cov[sapply(cov, function(x) length(x@values)) > 1] # in case some chr are empty, kick out those

### Initial ATAC-seq peak calling (IF peak file provided, skip)
if (is.character(peak_fn)) {
    ### user-defined peak regions, withou peak calling
    if (opt$verbose) {
        cat(">>> Read in input ATAC-seq peaks ...\n")
    }
    peaks <- tryCatch(
        import(peak_fn),
        error = function(e) {
            stop("Unable to read peak file.")
        }
    )
    ### fit to nb distribution
    mu <- sapply(cov, function(x) {
        cut <- quantile(x, 0.99) # kick out potential outliers
        x <- x[x <= cut]
        return(list(mu = mean(x), vi = var(x)))
    })
    vi <- unlist(mu[2, ])
    mu <- unlist(mu[1, ])
    vi <- ifelse(vi > mu, vi, mu + 1e-6)
    r <- mu * mu / (vi - mu)
} else if (is.null(peak_fn)) {
    ### call peaks using built-in peak caller
    if (opt$verbose) {
        cat(">>> Calling ATAC-seq peaks using RegBay's peak caller ...\n")
    }
    ### fit to nb distribution
    mu <- sapply(cov, function(x) {
        cut <- quantile(x, 0.99) # kick out potential outliers
        x <- x[x <= cut]
        return(list(mu = mean(x), vi = var(x)))
    })
    vi <- unlist(mu[2, ])
    mu <- unlist(mu[1, ])
    vi <- ifelse(vi > mu, vi, mu + 1e-6)
    r <- mu * mu / (vi - mu)
    peak.cut <- stats::qnbinom(1 - opt$pval, size = r, mu = mu, lower.tail = T)
    # call peaks using NB distribution
    cuts <- lapply(seq_along(cov), function(i) {
        IRanges::slice(cov[[i]], lower = peak.cut[i])
    })
    peaks <- mclapply(seq_along(cuts), function(i) {
        GRanges(seqnames = names(cov)[i], ranges = cuts[[i]]@ranges)
    }, mc.cores = opt$cores)
    peaks <- do.call("c", peaks)
    peaks <- peaks[width(peaks) > 100]
}


### cut long peaks (> 1kb) into individual peaks
if (opt$verbose) {
    cat(">>> Narrowing and spliting long peaks ...\n")
}
lens <- max(opt$cutoff, quantile(width(peaks), 0.8))
# shorts <- peaks[which(width(peaks) <= lens & width(peaks) > 100)]
cuts <- Views(cov, peaks)
longs <- lapply(cuts, function(x) x[width(x) > lens])

### Split-and-call for long peaks
long.slices <- lapply(longs, function(x) {
    slopes <- viewApply(x, diff)
    slices <- mclapply(seq_along(slopes), function(i) {
        ### Spliting - cut long peaks using derivatives
        y <- cumsum(slopes[[i]])
        y <- as.numeric(y)
        # smoothing
        y <- stats::filter(y, rep(1 / length(y) * 40, length(y) / 40), sides = 2, method = "convolution", circular = T)
        y <- diff(y)
        climbs <- which(y > 0.4)
        drops <- which(y < -0.4)
        # find the climbing segments
        ends <- which(diff(climbs) > 20)
        ends <- c(ends, length(climbs))
        climbs <- lapply(seq_along(ends), function(j) {
            if (j == 1) {
                vec <- climbs[1:ends[1]]
            } else {
                vec <- climbs[(ends[j - 1] + 1):ends[j]]
            }
            if (length(vec) < 20) {
                return()
            }
            return(vec)
        })
        climbs <- Filter(Negate(is.null), climbs)
        # find the dropping segments
        ends <- which(diff(drops) > 20)
        ends <- c(ends, length(drops))
        drops <- lapply(seq_along(ends), function(j) {
            if (j == 1) {
                vec <- drops[1:ends[1]]
            } else {
                vec <- drops[(ends[j - 1] + 1):ends[j]]
            }
            if (length(vec) < 20) {
                return()
            }
            return(vec)
        })
        drops <- Filter(Negate(is.null), drops)
        # merging to peaks
        slices <- list()
        if (length(climbs) > 0 & length(drops) > 0) {
            ends <- unlist(lapply(drops, function(x) x[1]))
            for (j in seq_along(climbs)) {
                k <- climbs[[j]]
                lefts <- k[length(k)]
                lefts <- c(lefts:(lefts + 80))
                pick <- unique(c(which(ends %in% lefts), which(sapply(drops, function(x) {
                    any(lefts %in% x)
                }))))
                if (length(pick) > 0) {
                    if (length(pick) == 1) {
                        vec <- c(k, lefts, drops[[pick]])
                    } else {
                        vec <- c(k, lefts, unlist(drops[pick]))
                    }
                    vec <- unique(vec)
                    slices <- append(slices, list(vec))
                } else {
                    next
                }
            }
            if (length(slices) == 0) {
                slices <- list(min(unlist(climbs)):max(unlist(drops)))
            }
        }

        ### peak calling based on 10 kb window
        segs <- x@ranges[i]
        start(segs) <- max(1, start(segs) - opt$window / 2)
        ind <- start(segs)
        end(segs) <- min(length(x@subject), end(segs) + opt$window / 2)
        segs <- Views(x@subject, segs)
        segs <- as.numeric(segs[[1]])
        v <- var(segs)
        k <- mean(segs)
        v <- ifelse(v < k, k, v)
        v <- k * k / (v - k + 1e-06)
        v <- 1 - stats::pnbinom(segs, size = v, mu = k, lower.tail = T)
        v <- v[(start(x@ranges[i]) - ind + 1):(end(x@ranges[i]) - ind + 1)]
        k <- which(v < opt$pval)
        if (length(slices) == 0) {
            v <- k
        } else {
            v <- k[-which(k %in% unlist(slices))] # remove the positions that are in slices
        }
        if (length(v) == 0) {
            return(slices)
        }
        ends <- which(diff(v) > 20)
        if (length(ends) == 0) {
            slices <- append(slices, list(min(v):max(v)))
        } else {
            ends <- c(ends, length(v))
            drops <- lapply(seq_along(ends), function(j) {
                if (j == 1) {
                    vec <- v[1:ends[1]]
                } else {
                    vec <- v[(ends[j - 1] + 1):ends[j]]
                }
                if (length(vec) < 20) {
                    return()
                }
                return(vec)
            })
            drops <- Filter(Negate(is.null), drops)
            slices <- append(slices, drops)
        }
        return(slices)
    }, mc.cores = opt$cores)
    filts <- sapply(slices, length)
    filts <- which(filts == 0)
    starts <- sapply(slices, function(y) sapply(y, min))
    ends <- sapply(slices, function(y) sapply(y, max))
    sub_views <- lapply(seq_along(x), function(i) {
        if (i %in% filts) {
            return(NULL)
        } else {
            return(IRanges(start = start(x[i]) + starts[[i]], end = start(x[i]) + ends[[i]]))
        }
    })
    sub_views <- Filter(function(x) length(x) > 0, sub_views)
    subs <- do.call("c", sub_views)
    subs <- reduce(subs)
    gaps <- gaps(subs)
    picks <- which(width(gaps) < 50)
    end(subs)[picks] <- end(subs)[picks + 1]
    subs <- reduce(subs)
    sub_views <- Views(subject(x), subs)
    return(sub_views)
})
longs <- mclapply(seq_along(long.slices), function(i) {
    GRanges(seqnames = names(cov)[i], ranges = long.slices[[i]]@ranges)
}, mc.cores = opt$cores)
longs <- do.call("c", longs)
longs <- longs[width(longs) > 100]


### local test for short peaks
cuts <- Views(cov, peaks)
shorts <- lapply(cuts, function(x) x[width(x) <= lens & width(x) > 100])
shorts <- lapply(shorts, function(x) { # skip the peaks at the first and last 5kb of each chromosome
    x[start(x) > 5000 | end(x) < length(x@subject) - 5000]
})
short.slices <- mclapply(shorts, function(x) {
    regs <- x@ranges
    start(regs) <- start(regs) - opt$window / 2
    end(regs) <- end(regs) + opt$window / 2

    deps <- Views(x@subject, regs)
    slices <- viewApply(deps, function(run) {
        vrun <- var(run)
        krun <- mean(run)
        vrun <- ifelse(vrun > krun, vrun, krun + 1e-6)
        vrun <- krun * krun / (vrun - krun)
        prun <- 1 - stats::pnbinom(as.numeric(run), size = vrun, mu = krun, lower.tail = T)
        prun <- prun[(opt$window / 2 + 1):(length(prun) - opt$window / 2)]
        krun <- which(prun < opt$pval)
        if (length(krun) == 0 | max(krun) - min(krun) < 100) {
            return(integer(0))
        } else {
            return(min(krun):max(krun))
        }
    })
    filts <- sapply(slices, length)
    filts <- which(filts == 0)
    starts <- sapply(slices, min)
    ends <- sapply(slices, max)
    sub_views <- lapply(seq_along(x), function(i) {
        if (i %in% filts) {
            return(NULL)
        } else {
            return(IRanges(start = start(x[i]) + starts[i], end = start(x[i]) + ends[i]))
        }
    })
    sub_views <- Filter(function(x) length(x) > 0, sub_views)
    subs <- do.call("c", sub_views)
    subs <- reduce(subs)
    gaps <- gaps(subs)
    picks <- which(width(gaps) < 50)
    end(subs)[picks] <- end(subs)[picks + 1]
    subs <- reduce(subs)
    sub_views <- Views(subject(x), subs)
    return(sub_views)
}, mc.cores = opt$cores)

lefts <- lapply(seq_along(shorts), function(i) {
    sigs <- short.slices[[i]]
    ori <- shorts[[i]]
    hits <- findOverlaps(sigs@ranges, ori@ranges, type = "any")
    return(ori[-unique(hits@to)])
})
left.slices <- lapply(lefts, function(x) {
    slopes <- viewApply(x, diff)
    slices <- mclapply(seq_along(slopes), function(i) {
        ### Spliting - cut long peaks using derivatives
        y <- cumsum(slopes[[i]])
        y <- as.numeric(y)
        # smoothing
        y <- stats::filter(y, rep(1 / length(y) * 40, length(y) / 40), sides = 2, method = "convolution", circular = T)
        y <- diff(y)
        climbs <- which(y > 0.4)
        drops <- which(y < -0.4)
        # find the climbing segments
        ends <- which(diff(climbs) > 20)
        ends <- c(ends, length(climbs))
        climbs <- lapply(seq_along(ends), function(j) {
            if (j == 1) {
                vec <- climbs[1:ends[1]]
            } else {
                vec <- climbs[(ends[j - 1] + 1):ends[j]]
            }
            if (length(vec) < 20) {
                return()
            }
            return(vec)
        })
        climbs <- Filter(Negate(is.null), climbs)
        # find the dropping segments
        ends <- which(diff(drops) > 20)
        ends <- c(ends, length(drops))
        drops <- lapply(seq_along(ends), function(j) {
            if (j == 1) {
                vec <- drops[1:ends[1]]
            } else {
                vec <- drops[(ends[j - 1] + 1):ends[j]]
            }
            if (length(vec) < 20) {
                return()
            }
            return(vec)
        })
        drops <- Filter(Negate(is.null), drops)
        # merging to peaks
        slices <- list()
        if (length(climbs) > 0 & length(drops) > 0) {
            ends <- unlist(lapply(drops, function(x) x[1]))
            for (j in seq_along(climbs)) {
                k <- climbs[[j]]
                lefts <- k[length(k)]
                lefts <- c(lefts:(lefts + 80))
                pick <- unique(c(which(ends %in% lefts), which(sapply(drops, function(x) {
                    any(lefts %in% x)
                }))))
                if (length(pick) > 0) {
                    vec <- c(k, lefts, drops[[pick]])
                    vec <- unique(vec)
                    slices <- append(slices, list(vec))
                } else {
                    next
                }
            }
            if (length(slices) == 0) {
                slices <- list(min(unlist(climbs)):max(unlist(drops)))
            }
        }
        return(slices)
    }, mc.cores = opt$cores)
    filts <- sapply(slices, length)
    filts <- which(filts == 0)
    starts <- sapply(slices, function(y) sapply(y, min))
    ends <- sapply(slices, function(y) sapply(y, max))
    sub_views <- lapply(seq_along(x), function(i) {
        if (i %in% filts) {
            return(NULL)
        } else {
            return(IRanges(start = start(x[i]) + starts[[i]], end = start(x[i]) + ends[[i]]))
        }
    })
    sub_views <- Filter(function(x) length(x) > 0, sub_views)
    if (length(sub_views) == 0) {
        return(NULL)
    }
    subs <- do.call("c", sub_views)
    subs <- reduce(subs)
    gaps <- gaps(subs)
    picks <- which(width(gaps) < 50)
    end(subs)[picks] <- end(subs)[picks + 1]
    subs <- reduce(subs)
    sub_views <- Views(subject(x), subs)
    return(sub_views)
})
lefts <- mclapply(seq_along(left.slices), function(i) {
    if (is.null(left.slices[[i]])) {
        return(NULL)
    }
    GRanges(seqnames = names(cov)[i], ranges = left.slices[[i]]@ranges)
}, mc.cores = opt$cores)
lefts <- do.call("c", lefts)
shorts <- mclapply(seq_along(short.slices), function(i) {
    GRanges(seqnames = names(cov)[i], ranges = short.slices[[i]]@ranges)
}, mc.cores = opt$cores)
shorts <- do.call("c", shorts)
shorts <- c(shorts, lefts)
shorts <- shorts[width(shorts) > 100]

peaks <- c(shorts, longs)
peaks <- reduce(peaks)
peaks <- split(peaks@ranges, peaks@seqnames)
peaks <- lapply(peaks, function(x) {
    gaps <- gaps(x)
    picks <- which(width(gaps) < 50)
    end(x[picks]) <- end(x[picks + 1])
    return(x)
})
peaks <- mclapply(seq_along(peaks), function(i) {
    GRanges(seqnames = names(cov)[i], ranges = peaks[[i]])
}, mc.cores = opt$cores)
peaks <- do.call("c", peaks)
peaks <- reduce(peaks)

### Calculate p- and q-values
raws <- Views(cov, peaks)
raws <- viewMaxs(raws)
pvalue <- lapply(raws, function(x) 1 - stats::pnbinom(x, size = r, mu = mu, lower.tail = T))
pvalue <- unlist(pvalue)
qvalue <- p.adjust(pvalue, method = "BH")
peaks$score <- qvalue
peaks <- peaks[peaks$score < opt$pval]

### output peak regions
export.bed(peaks, paste0(opt$outdir, '/peaks.bed'))
cat(paste0("Peak processing done! Run time: ", round(difftime(Sys.time(), s.time, units = "mins"), digits = 2), " min\n\n"))
