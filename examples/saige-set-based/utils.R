# utils.R
# Utility functions for SAIGE result visualization

#' Extract variant encoding parameters from VCF
#'
#' Reads the VCF file to extract genotype frequencies (r, h, a) and
#' rescaling parameters (globalMinDomDosage, globalMaxDomDosage) for a specific variant.
#'
#' @param vcf_file Path to the non-additive VCF file (gzipped)
#' @param variant_id Variant ID to extract parameters for
#' @return List with r, h, a, global_min, global_max
extract_vcf_params <- function(vcf_file, variant_id) {
    # Read the VCF line for this variant
    cmd <- sprintf("gunzip -c %s | grep -v '^##' | grep '%s'", vcf_file, variant_id)
    vcf_line <- system(cmd, intern = TRUE)

    if (length(vcf_line) == 0) {
        stop(sprintf("Variant %s not found in VCF file %s", variant_id, vcf_file))
    }

    # Extract INFO field (8th column)
    info_field <- strsplit(vcf_line, "\t")[[1]][8]

    # Parse INFO field to extract r, h, a, globalMinDomDosage, globalMaxDomDosage
    info_parts <- strsplit(info_field, ";")[[1]]

    get_info_value <- function(key) {
        val_str <- grep(paste0("^", key, "="), info_parts, value = TRUE)
        if (length(val_str) == 0) return(NA)
        as.numeric(sub(paste0("^", key, "="), "", val_str))
    }

    list(
        r = get_info_value("r"),
        h = get_info_value("h"),
        a = get_info_value("a"),
        global_min = get_info_value("globalMinDomDosage"),
        global_max = get_info_value("globalMaxDomDosage")
    )
}

#' Calculate dose-response relationship from SAIGE results
#'
#' This function implements the orthogonalized Gram-Schmidt encoding approach
#' to calculate genetic effects for each genotype dosage (0, 1, 2). The encoding
#' uses genotype frequencies from the VCF and applies the same linear rescaling
#' that was used before SAIGE fitting.
#'
#' The calculation follows these steps:
#'   1. Extract genotype frequencies (r, h, a) from VCF
#'   2. Compute orthogonalized X^A and X^D using Gram-Schmidt process
#'   3. Apply linear rescaling to X^D using global min/max from VCF
#'   4. Calculate effects relative to dosage 0 (reference)
#'
#' @param beta_add Additive effect estimate (from SAIGE additive model)
#' @param se_add Standard error of additive effect
#' @param beta_dom Dominance deviation effect estimate (from SAIGE non-additive model)
#' @param se_dom Standard error of dominance deviation effect
#' @param vcf_file Path to the non-additive VCF file used for SAIGE
#' @param variant_id Variant ID to calculate dose-response for
#' @return Data frame with columns: Dosage, X_A, X_D, Estimate, SE, Lower_CI, Upper_CI
calculate_dose_response <- function(beta_add, se_add, beta_dom, se_dom,
                                   vcf_file, variant_id) {

    # Extract parameters from VCF
    params <- extract_vcf_params(vcf_file, variant_id)
    r <- params$r
    h <- params$h
    a <- params$a
    min_dosage <- params$global_min
    max_dosage <- params$global_max

    # Calculate orthogonalized additive encoding (Equation 5 from paper)
    # X^A = [-(h + 2a), 1-(h + 2a), 2-(h + 2a)] / sqrt(a+r - (a-r)^2)
    denom_a <- sqrt(a + r - (a - r)^2)
    x_a <- c(-(h + 2*a), 1 - (h + 2*a), 2 - (h + 2*a)) / denom_a

    # Calculate orthogonalized dominance deviation encoding (Equation 6 from paper)
    # X^D = [-ha, 2ar, -hr] / sqrt(har(a+r - (a-r)^2))
    # These are the ORIGINAL values before rescaling
    scaling_non_add <- 1 / sqrt(h * a * r * (a + r - (a - r)^2))
    x_d_orig <- scaling_non_add * c(-h*a, 2*a*r, -h*r)

    # Rescale beta_dom back to the original scale (before VCF rescaling)
    # SAIGE fits on rescaled X^D, so we need to transform beta back
    beta_dom_rescaled <- beta_dom * (2 / ((max_dosage - min_dosage) * scaling_non_add))
    se_dom_rescaled <- se_dom * (2 / ((max_dosage - min_dosage) * scaling_non_add))

    # Calculate absolute effects for each dosage using ORIGINAL X^D values
    # Equation: g_hat = beta_a * X^A + beta_d_original * X^D_original
    effects_abs <- beta_add * x_a + beta_dom_rescaled * x_d_orig

    # Make effects relative to dosage 0 (reference genotype)
    effects <- effects_abs - effects_abs[1]

    # Calculate standard errors relative to dosage 0
    # SE(g_d - g_0) = sqrt((X^A_d - X^A_0)^2 * SE_a^2 + (X^D_d - X^D_0)^2 * SE_d_rescaled^2)
    x_a_diff <- x_a - x_a[1]
    x_d_diff <- x_d_orig - x_d_orig[1]
    ses <- sqrt(x_a_diff^2 * se_add^2 + x_d_diff^2 * se_dom_rescaled^2)

    # Calculate 99% confidence intervals
    lower_ci <- effects - 2.576 * ses
    upper_ci <- effects + 2.576 * ses

    # Return results as data frame
    data.frame(
        Dosage = 0:2,
        X_A = x_a,
        X_D = x_d_orig,  # Return original (unrescaled) X^D values
        Estimate = effects,
        SE = ses,
        Lower_CI = lower_ci,
        Upper_CI = upper_ci
    )
}

#' Plot dose-response relationship
#'
#' Creates a visualization showing the estimated genetic effect for each
#' genotype dosage (0, 1, 2) with 99% confidence intervals.
#'
#' @param res Results data frame from calculate_dose_response
#' @param variant_id Variant identifier for plot title
#' @param p_add P-value from additive test
#' @param p_dom P-value from dominance test
#' @param annotation Optional annotation label (e.g., "pLoF", "synonymous")
plot_dose_response <- function(res, variant_id, p_add, p_dom, annotation = NULL) {

    # Setup plot margins and ranges
    y_min <- min(res$Lower_CI)
    y_max <- max(res$Upper_CI)
    range <- y_max - y_min
    y_lim <- c(y_min - 0.15 * range, y_max + 0.15 * range)

    # Create base plot
    plot(res$Dosage, res$Estimate,
         type = "n",
         ylim = y_lim,
         xlim = c(-0.2, 2.2),
         xaxt = "n",
         xlab = "Genotype Dosage",
         ylab = "Estimated Effect",
         main = variant_id)

    # Add subtitle with P-values and annotation
    subtitle <- sprintf("P_add: %.2e  |  P_dom: %.2e", p_add, p_dom)
    if (!is.null(annotation)) {
        subtitle <- paste0(subtitle, "  |  ", annotation)
    }
    mtext(subtitle, side = 3, line = 0.2, cex = 0.8)

    # Add grid and reference line
    grid(nx = NA, ny = NULL, col = "gray90")
    abline(h = 0, col = "black", lty = 2, lwd = 1)

    # Custom x-axis labels
    axis(1, at = 0:2, labels = c("Hom Ref (0)", "Het (1)", "Hom Alt (2)"))

    # Plot 99% CI as shaded region
    polygon(c(res$Dosage, rev(res$Dosage)),
            c(res$Lower_CI, rev(res$Upper_CI)),
            col = rgb(0.2, 0.4, 0.8, 0.2),
            border = NA)

    # Plot effect estimates
    lines(res$Dosage, res$Estimate, col = "firebrick", lwd = 2)
    points(res$Dosage, res$Estimate, pch = 19, col = "firebrick", cex = 1.5)

    # Add legend
    legend("topleft",
           legend = c("Estimate", "99% CI"),
           col = c("firebrick", rgb(0.2, 0.4, 0.8, 0.2)),
           lwd = c(2, 10),
           pch = c(19, NA),
           bty = "n")
}

#' Plot dose-response comparison between two variants
#'
#' Creates a single plot comparing dose-response curves for two variants,
#' typically used to contrast causal (e.g., pLoF) vs null (e.g., synonymous) variants.
#'
#' @param d1 Data row for first variant (from merged SAIGE results)
#' @param d2 Data row for second variant
#' @param vcf_file Path to non-additive VCF file
#' @param label1 Label for first variant (e.g., "pLoF")
#' @param label2 Label for second variant (e.g., "Synonymous")
#' @param color1 Color for first variant (default: "#d62728" red)
#' @param color2 Color for second variant (default: "#2ca02c" green)
plot_dose_response_comparison <- function(d1, d2, vcf_file,
                                         label1 = "Variant 1", label2 = "Variant 2",
                                         color1 = "#d62728", color2 = "#2ca02c") {

    # Calculate dose-response for both variants
    dose1 <- calculate_dose_response(
        beta_add = d1$BETA.add, se_add = d1$SE.add,
        beta_dom = d1$BETA.dom, se_dom = d1$SE.dom,
        vcf_file = vcf_file, variant_id = d1$MarkerID
    )

    dose2 <- calculate_dose_response(
        beta_add = d2$BETA.add, se_add = d2$SE.add,
        beta_dom = d2$BETA.dom, se_dom = d2$SE.dom,
        vcf_file = vcf_file, variant_id = d2$MarkerID
    )

    # Setup plot range
    y_range <- range(c(dose1$Lower_CI, dose1$Upper_CI,
                      dose2$Lower_CI, dose2$Upper_CI))
    y_pad <- diff(y_range) * 0.15
    y_lim <- y_range + c(-y_pad, y_pad)

    # Create base plot
    plot(NULL, xlim = c(-0.2, 2.2), ylim = y_lim, xaxt = "n",
         xlab = "Genotype Dosage", ylab = "Estimated Effect",
         main = sprintf("Dose-Response: %s vs %s", label1, label2))

    # Grid and reference
    grid(nx = NA, ny = NULL, col = "gray90")
    abline(h = 0, col = "black", lty = 2, lwd = 1)
    axis(1, at = 0:2, labels = c("Hom Ref (0)", "Het (1)", "Hom Alt (2)"))

    # Plot first variant
    col1_rgb <- col2rgb(color1) / 255
    polygon(c(dose1$Dosage, rev(dose1$Dosage)),
            c(dose1$Lower_CI, rev(dose1$Upper_CI)),
            col = rgb(col1_rgb[1], col1_rgb[2], col1_rgb[3], 0.2), border = NA)
    lines(dose1$Dosage, dose1$Estimate, col = color1, lwd = 2.5)
    points(dose1$Dosage, dose1$Estimate, pch = 19, col = color1, cex = 1.5)

    # Plot second variant
    col2_rgb <- col2rgb(color2) / 255
    polygon(c(dose2$Dosage, rev(dose2$Dosage)),
            c(dose2$Lower_CI, rev(dose2$Upper_CI)),
            col = rgb(col2_rgb[1], col2_rgb[2], col2_rgb[3], 0.2), border = NA)
    lines(dose2$Dosage, dose2$Estimate, col = color2, lwd = 2.5)
    points(dose2$Dosage, dose2$Estimate, pch = 19, col = color2, cex = 1.5)

    # Legend
    legend("topleft",
           legend = c(sprintf("%s (P=%.1e)", label1, d1$p.value.add),
                     sprintf("%s (P=%.1e)", label2, d2$p.value.add),
                     "99% CI"),
           col = c(color1, color2, "gray"),
           lwd = c(2.5, 2.5, 10),
           pch = c(19, 19, NA),
           bty = "n")
}

#' Get expected P-values for QQ plot
#'
#' @param p_values Vector of observed P-values
#' @return Vector of expected P-values under uniform distribution
get_expected_p <- function(p_values) {
    n <- length(p_values)
    rank <- rank(p_values, ties.method = "first")
    expected <- (rank - 0.5) / n
    return(expected)
}
