# utils.R
# Utility functions for GWAS analysis

#' Calculate dose-response relationship
#'
#' @param model_add Fitted additive model (lm object)
#' @param model_dom Fitted dominance model (lm object)
#' @param data Data table containing IID, Additive, Dominance columns
#' @param variant_id ID of the variant to analyze
#' @return Data frame with estimates and CI for dosages 0, 1, 2
calculate_dose_response <- function(model_add, model_dom, data, variant_id) {
    # Extract coefficients and SEs
    # Note: We assume the models are: Y ~ Additive and Y ~ Dominance
    
    coef_add <- summary(model_add)$coefficients
    beta_a <- coef_add["Additive", "Estimate"]
    se_a <- coef_add["Additive", "Std. Error"]
    
    coef_dom <- summary(model_dom)$coefficients
    # It's possible the dominance model failed or has no term if constant
    if ("Dominance" %in% rownames(coef_dom)) {
        beta_d <- coef_dom["Dominance", "Estimate"]
        se_d <- coef_dom["Dominance", "Std. Error"]
    } else {
        beta_d <- 0
        se_d <- 0
    }
    
    # Identify X^D values for each additive dosage d in {0, 1, 2}
    # We infer this from the data.
    # We act on the specific variant subset passed in 'data'
    
    dosages <- c(0, 1, 2)
    results <- data.frame(
        Dosage = dosages,
        X_A = dosages, # Additive genotype is definitionally 0, 1, 2
        X_D = NA,      # To be determined
        Estimate = NA,
        SE = NA,
        Lower_CI = NA,
        Upper_CI = NA
    )
    
    for (i in 1:length(dosages)) {
        d <- dosages[i]
        
        # Find samples with this additive dosage (approximate to handle float)
        # We allow small tolerance
        subset_mask <- abs(data$Additive - d) < 0.1
        
        if (sum(subset_mask) > 0) {
            # Get mean Dominance dosage for these samples
            x_d_val <- mean(data$Dominance[subset_mask], na.rm=TRUE)
        } else {
            # If no samples with this genotype, we can't empirically determine X^D
            # For visualization, we might mark as NA
            x_d_val <- NA
        }
        
        results$X_D[i] <- x_d_val
        
        if (!is.na(x_d_val)) {
            # Equation 1: g_hat = beta_a * X_a + beta_d * X_d
            g_hat <- beta_a * d + beta_d * x_d_val
            
            # Equation 3 (implied): SE = sqrt((X_a * SE_a)^2 + (X_d * SE_d)^2)
            # Covariance is assumed 0 due to orthogonality
            se_full <- sqrt((d * se_a)^2 + (x_d_val * se_d)^2)
            
            results$Estimate[i] <- g_hat
            results$SE[i] <- se_full
            results$Lower_CI[i] <- g_hat - 1.96 * se_full
            results$Upper_CI[i] <- g_hat + 1.96 * se_full
        }
    }
    
    return(results)
}

#' Print dose-response table
print_dose_response <- function(res, variant_id) {
    cat(sprintf("\nDose-Response for Variant: %s\n", variant_id))
    cat("Dosage (d) | X^A |   X^D   | Estimate (g_hat) |    SE    |   95% CI   \n")
    cat("-----------+-----+---------+------------------+----------+------------\n")
    for (i in 1:nrow(res)) {
        if (!is.na(res$Estimate[i])) {
            cat(sprintf("     %d     |  %d  | %7.4f |     %7.4f      | %8.4f | [%.3f, %.3f]\n",
                        res$Dosage[i], res$X_A[i], res$X_D[i], 
                        res$Estimate[i], res$SE[i], 
                        res$Lower_CI[i], res$Upper_CI[i]))
        } else {
            cat(sprintf("     %d     |  %d  |    NA   |        NA        |    NA    |      NA     \n",
                        res$Dosage[i], res$X_A[i]))
        }
    }
}

#' Plot dose-response relationship
#'
#' @param res Results data frame from calculate_dose_response
#' @param variant_id Title for the plot
#' @param phenotype Phenotype name
#' @param p_add Additive P-value
#' @param p_dom Dominance P-value
plot_dose_response <- function(res, variant_id, phenotype, p_add, p_dom) {
    # Remove NAs for plotting
    df <- res[!is.na(res$Estimate), ]
    if (nrow(df) < 2) return()
    
    # Setup plot range
    y_min <- min(df$Lower_CI)
    y_max <- max(df$Upper_CI)
    range <- y_max - y_min
    # Add some padding
    y_lim <- c(y_min - 0.1 * range, y_max + 0.1 * range)
    
    # Base Plot
    plot(df$Dosage, df$Estimate, 
         type="n", # No points yet
         ylim=y_lim, 
         xlim=c(-0.2, 2.2),
         xaxt="n",
         xlab="Genotype Dosage (d)", 
         ylab="Estimated Effect (g_hat)",
         main=paste0(phenotype, "\nVariant: ", variant_id))
    
    # Add P-values as Subtitle
    mtext(sprintf("P_add: %.2e  |  P_dom: %.2e", p_add, p_dom), side=3, line=0.2, cex=0.8)
    
    # Grid
    grid(nx=NA, ny=NULL, col="gray90")
    
    # Reference Line at Y=0
    abline(h=0, col="black", lty=2)
    
    # X Axis Labels
    axis(1, at=0:2, labels=c("Hom Ref (0)", "Het (1)", "Hom Alt (2)"))
    
    # Shaded 95% CI
    polygon(c(df$Dosage, rev(df$Dosage)), 
            c(df$Lower_CI, rev(df$Upper_CI)),
            col=rgb(0.2, 0.4, 0.8, 0.2), # Transparent Blue
            border=NA)
            
    # Lines and Points
    lines(df$Dosage, df$Estimate, col="firebrick", lwd=2, lty=1)
    points(df$Dosage, df$Estimate, pch=19, col="firebrick", cex=1.5)
    
    # Add Legend
    legend("topleft", legend=c("Estimate", "95% CI"), 
           col=c("firebrick", rgb(0.2, 0.4, 0.8, 0.2)), 
           lwd=c(2, 10), pch=c(19, NA), bty="n")
}
