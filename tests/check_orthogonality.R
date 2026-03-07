#!/usr/bin/env Rscript

# Script to verify orthogonality between Additive and Dominance VCFs
# Usage: Rscript check_orthogonality.R --additive <vcf> --dominance <vcf>

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript check_orthogonality.R --additive <vcf> --dominance <vcf>")
}

add_file <- args[which(args == "--additive") + 1]
dom_file <- args[which(args == "--dominance") + 1]

if (!file.exists(add_file)) stop(paste("Additive file not found:", add_file))
if (!file.exists(dom_file)) stop(paste("Dominance file not found:", dom_file))

# Function to read VCF dosages
read_dosages <- function(file, is_dominance=FALSE) {
    # Use grep to skip header
    cmd <- paste("gunzip -c", file, "| grep -v '^#'")
    data <- tryCatch(read.table(text=system(cmd, intern=TRUE), stringsAsFactors=FALSE), error=function(e) NULL)
    
    if (is.null(data) || nrow(data) == 0) {
        return(list())
    }
    
    # Cols 1-9 are metadata. 10+ are samples.
    # Format col is 9.
    
    res <- list()
    for (i in 1:nrow(data)) {
        var_id <- data[i, 3] # ID column
        format_str <- data[i, 9]
        raw_vals <- data[i, 10:ncol(data)]
        
        # Parse dosages
        # If Dominance: Look for DS
        # If Additive: Look for GT and convert to 0,1,2
        
        format_parts <- strsplit(format_str, ":")[[1]]
        
        if (is_dominance) {
            ds_idx <- which(format_parts == "DS")
            if (length(ds_idx) > 0) {
                 vals <- as.numeric(sapply(strsplit(as.character(raw_vals), ":"), function(x) x[ds_idx]))
            } else {
                 # Fallback if just raw numbers (unlikely for VCF but possible in some formats)
                 vals <- as.numeric(raw_vals)
            }
        } else {
            # Additive - assume GT
            gt_idx <- which(format_parts == "GT")
            if (length(gt_idx) > 0) {
                vals <- sapply(strsplit(as.character(raw_vals), ":"), function(x) x[gt_idx])
                # Convert 0|0 -> 0, 0|1 -> 1, 1|1 -> 2
                vals <- vapply(vals, function(x) {
                    if (grepl("\\.", x)) return(NA_real_)
                    parts <- strsplit(x, "[|/]")[[1]]
                    sum(as.numeric(parts))
                }, numeric(1))
            } else {
                # Try DS if exists
                ds_idx <- which(format_parts == "DS")
                if (length(ds_idx) > 0) {
                    vals <- as.numeric(sapply(strsplit(as.character(raw_vals), ":"), function(x) x[ds_idx]))
                } else {
                    stop(paste("Cannot parse additive dosage for variant", var_id))
                }
            }
        }
        res[[var_id]] <- vals
    }
    return(res)
}

cat("Reading Additive VCF...\n")
add_data <- read_dosages(add_file, is_dominance=FALSE)

cat("Reading Dominance VCF...\n")
dom_data <- read_dosages(dom_file, is_dominance=TRUE)

# Check Variants
common_vars <- intersect(names(add_data), names(dom_data))
cat(sprintf("Found %d common variants.\n", length(common_vars)))

epsilon <- 1e-4

for (v in common_vars) {
    A <- add_data[[v]]
    D <- dom_data[[v]]
    
    # Remove missing
    valid <- !is.na(A) & !is.na(D)
    A <- A[valid]
    D <- D[valid]
    
    if (length(A) == 0) next
    
    # Calculate Centered Dot Product (Covariance-like check)
    # Statistical orthogonality means Cov(A, D) ~ 0
    # This is equivalent to cor(A, D) ~ 0 or sum((A - mean(A)) * (D - mean(D))) ~ 0
    
    # We use correlation if variance is non-zero, else 0
    sd_a <- sd(A)
    sd_d <- sd(D)
    
    if (sd_a == 0 || sd_d == 0) {
        val <- 0 # Monomorphic is orthogonal to anything by constant definition
    } else {
        # Check sum of centered product to be close to 0
        # cor() might be NaN if extremely constant, but sd check handles it
        val <- cor(A, D)
    }
    
    val[is.na(val)] <- 0
    
    cat(sprintf("Variant %s: Correlation = %.6f\n", v, val))
    
    if (abs(val) > epsilon) {
        stop(sprintf("FAILURE: Variant %s is NOT statistically orthogonal. Correlation = %f", v, val))
    }
}

cat("\nSUCCESS: All variants are statistically orthogonal (uncorrelated).\n")
