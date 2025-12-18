# Determine script directory to load utils
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)])
if (length(script_name) == 0) { 
  script_dir <- "examples" # Fallback
} else {
  script_dir <- dirname(script_name)
}

# Load utils
utils_path <- file.path(script_dir, "utils.R")
if (!file.exists(utils_path)) {
    # Try looking in examples/ relative to CWD if script_dir was just "." or similar logic failed
    if (file.exists("examples/utils.R")) {
        utils_path <- "examples/utils.R"
    } else {
        stop("Could not find utils.R")
    }
}
source(utils_path)

# Parse Args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript run_gwas.R --additive <vcf> --dominance <vcf> --phenotype <file> --out <prefix>")
}

add_file <- args[which(args == "--additive") + 1]
dom_file <- args[which(args == "--dominance") + 1]
pheno_file <- args[which(args == "--phenotype") + 1]
out_prefix <- args[which(args == "--out") + 1]

# Load Packages
if (!require("data.table")) install.packages("data.table", repos="https://cloud.r-project.org")
library(data.table)

# Function to read VCF dosage
read_vcf_dosage <- function(vcf_file) {
    # Check if file exists
    if (!file.exists(vcf_file)) stop(paste("File not found:", vcf_file))
    
    # Read header
    cmd_header <- paste("gunzip -c", vcf_file, "| grep '#CHROM' | head -n 1")
    header_line <- system(cmd_header, intern=TRUE)
    if (length(header_line) == 0) stop("Could not read VCF header")
    header <- strsplit(header_line, "\t")[[1]]
    samples <- header[10:length(header)]
    
    # Check number of lines (variants)
    # For now, we assume we want to analyze ALL variants in the file.
    # But if there are many, we should loop. 
    # For the example, we'll read the first variant or assume small file.
    # The example creates 2 variants or 1 gene.
    
    cmd <- paste("gunzip -c", vcf_file, "| grep -v '##' | grep -v '#CHROM'")
    data <- fread(cmd = cmd, header=FALSE)
    if (nrow(data) == 0) return(NULL)
        
    # We will reshape this to long format for analysis: IID, Variant, Dosage
    # But for the simplest example, let's just analyze the FIRST variant or GENE present.
    # If multiple variants, we run analysis for each.
    
    res_list <- list()
    
    for (i in 1:nrow(data)) {
        # Extract ID (Column 3)
        var_id <- as.character(data[i, 3])
        
        # Parse FORMAT (Column 9)
        format_str <- as.character(data[i, 9])
        format_parts <- strsplit(format_str, ":")[[1]]
        ds_idx <- which(format_parts == "DS")
        
        # Get raw sample columns (10 to end)
        # unlist to get character vector
        raw_values <- unlist(data[i, 10:ncol(data), with=FALSE])
        
        if (length(ds_idx) > 0) {
            # DS is present
            if (format_str == "DS") {
                dosages <- as.numeric(raw_values)
            } else {
                dosages <- as.numeric(sapply(strsplit(raw_values, ":"), function(x) x[ds_idx]))
            }
        } else {
            # Check for GT
            gt_idx <- which(format_parts == "GT")
            if (length(gt_idx) > 0) {
                 # Extract GT part
                 if (format_str == "GT") {
                     gts <- raw_values
                 } else {
                     gts <- sapply(strsplit(raw_values, ":"), function(x) x[gt_idx])
                 }
                 # Calc dosage from GT (count 1s)
                 # e.g. 0|0 -> 0, 0|1 -> 1, 1|1 -> 2, .|. -> NA
                 # Simple approach: sum of integers? No, 0|1 is string.
                 # Robust way: count '1's.
                 dosages <- vapply(gts, function(x) {
                     if (grepl("\\.", x)) return(NA_real_)
                     # sum of numbers? 
                     # split by / or |
                     alleles <- strsplit(x, "[|/]")[[1]]
                     sum(as.numeric(alleles))
                 }, numeric(1))
            } else {
                 # Fallback: try direct numeric
                 dosages <- as.numeric(raw_values)
            }
        }
        
        res_list[[var_id]] <- data.table(IID=samples, Variant=var_id, Dosage=dosages)
    }
    
    return(rbindlist(res_list))
}

# 1. Load Phenotypes
pheno <- fread(pheno_file)

# 2. Load Genotypes
cat(sprintf("Loading Additive VCF: %s\n", add_file))
add_long <- read_vcf_dosage(add_file)
setnames(add_long, "Dosage", "Additive")

cat(sprintf("Loading Dominance VCF: %s\n", dom_file))
dom_long <- read_vcf_dosage(dom_file)
setnames(dom_long, "Dosage", "Dominance")

# Merge Genotypes
genotypes <- merge(add_long, dom_long, by=c("IID", "Variant"))

# Merge Phenotypes (All phenotypes: Y_rec, Y_add, Y_null)
df_all <- merge(genotypes, pheno, by="IID")

# Analysis Function
run_analysis <- function(df, y_col, title) {
    cat(sprintf("\n--- Phenotype: %s ---\n", y_col))
    
    # Loop over variants/genes
    variants <- unique(df$Variant)
    
    for (v in variants) {
        sub_df <- df[Variant == v]
        
        # Fit separate models for marginal estimates
        f_add <- as.formula(paste(y_col, "~ Additive"))
        f_dom <- as.formula(paste(y_col, "~ Dominance"))
        f_jnt <- as.formula(paste(y_col, "~ Additive + Dominance"))
        
        model_add <- lm(f_add, data=sub_df)
        model_dom <- lm(f_dom, data=sub_df)
        model_joint <- lm(f_jnt, data=sub_df)
        
        s <- summary(model_joint)
        coefs <- s$coefficients
        
        # Check if Dominance exists
        if ("Dominance" %in% rownames(coefs)) {
             est_dom <- coefs["Dominance", "Estimate"]
             p_dom <- coefs["Dominance", "Pr(>|t|)"]
        } else {
             est_dom <- NA
             p_dom <- NA
        }
        
        if ("Additive" %in% rownames(coefs)) {
             est_add <- coefs["Additive", "Estimate"]
             p_add <- coefs["Additive", "Pr(>|t|)"]
        } else {
             est_add <- NA
             p_add <- NA
        }
        
        cat(sprintf("Variant: %-10s | Additive: Est=%.2f P=%.2e | Dominance: Est=%.2f P=%.2e | R2=%.4f\n", 
                    v, est_add, p_add, est_dom, p_dom, s$r.squared))
                    
        # Calculate dose-response (only if decent size to avoid errors)
        if (nrow(sub_df) > 10) {
             dr <- calculate_dose_response(model_add, model_dom, sub_df, v)
             print_dose_response(dr, v)
             
             # Plot
             plot_dose_response(dr, v, title, p_add, p_dom)
        }
    }
}

# Run for all phenotypes
cat("\n=======================================================\n")
cat("  GWAS Analysis Results\n")
cat("=======================================================\n")

# Open PDF device
pdf_file <- paste0(out_prefix, "_plots.pdf")
pdf(pdf_file, width=8, height=6)
cat(sprintf("\nSaving plots to: %s\n", pdf_file))

for (y in c("Y_rec", "Y_add", "Y_null")) {
    run_analysis(df_all, y, y)
}

# Close device
dev.off()

cat("\nDone.\n")
