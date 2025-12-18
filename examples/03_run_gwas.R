# R Script to analyze generated data and demonstrate orthogonality

# Install/Load required packages
if (!require("data.table")) install.packages("data.table", repos="https://cloud.r-project.org")
library(data.table)

# Function to read VCF dosage (simplified)
read_vcf_dosage <- function(vcf_file) {
    # Using specific command to extract dosage/GT from our pseudo VCF
    # The MakePseudoVCF tool places dosages in the GT field or DS field depending on mode?
    # Let's inspect the output format in the bash script first or assume standard dosage extraction.
    # For make_pseudo_vcf, it usually outputs GTs like 0/0, 0/1, 1/1 which map to 0,1,2.
    
    cmd <- paste("gunzip -c", vcf_file, "| grep -v '##' | grep -v '#CHROM'")
    data <- fread(cmd = cmd, header=FALSE)
    
    # Header
    cmd_header <- paste("gunzip -c", vcf_file, "| grep '#CHROM'")
    header <- strsplit(system(cmd_header, intern=TRUE), "\t")[[1]]
    colnames(data) <- header
    
    # Convert GT to dosage
    # 0/0 -> 0, 0/1 -> 1, 1/1 -> 2
    # Our tools might output 1/0 or 0/1.
    
    # We want to return a matrix of Sample x Gene
    # Transpose
    samples <- header[10:length(header)]
    
    # For this simple example, we only have one gene 'GeneA'
    # Let's extract row for GeneA
    gene_row <- data[ID == "GeneA"]
    if (nrow(gene_row) == 0) return(NULL)
        
    # The columns 10:ncol(data) contain the dosages directly (e.g. 0, 1, 2 or floats)
    # We just need to unlist and convert to numeric
    dosages <- as.numeric(gene_row[, 10:ncol(data), with=FALSE])
    
    return(data.table(IID=samples, Dosage=dosages))
}

# 1. Load Phenotypes
pheno <- fread("examples/phenotypes.txt")

# 2. Load Genotypes
# Additive Encoding
add_dt <- read_vcf_dosage("examples/output/encoded_additive.vcf.gz")
setnames(add_dt, "Dosage", "Additive")

# Dominance (Non-additive) Encoding
dom_dt <- read_vcf_dosage("examples/output/encoded_dominance.vcf.gz")
setnames(dom_dt, "Dosage", "Dominance") 
# Note: In 'dominance' mode, the tool outputs a transformation where heterozygotes are emphasized 
# or residuals are captured strictly orthogonal to additive.
# Standard 0,1,2 additive coding.
# Orthogonal dominance often looks like: 0 -> 0, 1 -> 1, 2 -> 0 (if strictly het) or similar depending on implementation.
# Let's see what the regression says.

# Merge
df <- merge(pheno, add_dt, by="IID")
df <- merge(df, dom_dt, by="IID")

print("Data Summary:")
print(summary(df))

# 3. Correlation Check (Orthogonality)
cor_val <- cor(df$Additive, df$Dominance)
print(paste("Correlation between Additive and Dominance encodings:", round(cor_val, 4)))
# Ideally this should be low if orthogonalized, or at least they capture different information.

# 4. Regression Analysis
# Model: Y ~ Additive + Dominance
model <- lm(Y ~ Additive + Dominance, data=df)

print("Regression Results:")
print(summary(model))

# Interpretation
# If our phenotype is truly Recessive (effect only on 2 alleles), 
# we expect the 'Dominance' term to catch the deviation from the slope 
# that 'Additive' tries to fit.
# A purely additive model fits a straight line.
# A recessive model is 0, 0, 1.
# Additive sees: 0, 0.5, 1 (slope).
# Dominance sees: 0, 1, 0 (heterozygote) -> Recessive is -0.5 deviations? 
# OR if the tool outputs 0, 1, 0 for hets, then recessive (0,0,1) is 'low dominance'?
# Wait, let's verify what 'Dominance' encoding means in this tool.
# Usually:
# Genotype: 0, 1, 2
# Additive: 0, 1, 2
# Dominance: 0, 1, 0 (Classic)
# Orthogonal Dominance: 0, 2pq, -2p^2 (Genotyped based orthogonality)

# Let's check coefficients.
