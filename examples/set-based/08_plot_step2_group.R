
library(ggplot2)
library(data.table)
library(bravastring) # importing get_expected_p() from devtools::install_github("frhl/bravastring")
setwd("~/Projects/06_call_chets/call_chets/examples/set-based/")

# for pretty colors
col_palette <- c(
  "synonymous" = "#2ca02c",
  "pLoF" = "#d62728"
)

# gene-level results (note, because of a bug in SAIGE, we need to
# run each group separately for the non-additive encoding. You can
# run them together, but then the "Group" label will be messed up)
d1 <- fread("output/saige.step2.nonadditive.group.synonymous.txt")
d1 <- d1[d1$Group != "Cauchy"]
d2 <- fread("output/saige.step2.nonadditive.group.pLoF.txt")
d2 <- d2[d2$Group != "Cauchy"]
d <- rbind(d1, d2)

# keep only relevant cols
cols_to_keep <- c("Region", "Group", "Pvalue_Burden", "BETA_Burden", "SE_Burden", "Number_rare")
d <- d[,..cols_to_keep]

# merge with additive results (plof and synonymous)
da <- fread("output/saige.step2.additive.group.txt")
da <- da[Group != "Cauchy",..cols_to_keep]
d <- merge(d, da, all.x=TRUE, by=cols_to_keep[1:2], suffixes = c(".dom", ".add"))

# get expected P
d[, Pvalue_Burden.dom.expt := get_expected_p(Pvalue_Burden.dom)]
d[, Pvalue_Burden.add.expt := get_expected_p(Pvalue_Burden.add)]

# check pLoF
d[Group=="pLoF",]

# We also add a 'model' column to distinguish between them
plot_data <- rbind(
  d[, .(p.value = Pvalue_Burden.add, expt.p = Pvalue_Burden.add.expt, annotation = Group, model = "Additive [0,1,2]")],
  d[, .(p.value = Pvalue_Burden.dom, expt.p = Pvalue_Burden.dom.expt, annotation = Group, model = "Non-additive")]
)

# 2. Plot side-by-side
ggplot(plot_data, aes(x = -log10(expt.p), y = -log10(p.value), color = annotation)) +
  geom_point(data=plot_data[annotation != "pLoF"], size=2.5) +
  geom_point(data=plot_data[annotation == "pLoF"], size=2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values=col_palette) +
  facet_wrap(~model) + 
  theme_bw() +
  labs(
    x = "-log10(Expected P-value)",
    y = "-log10(Observed P-value)",
    color = "Annotation"
  ) + 
  ggtitle("Simulating a completely recessive trait", "Group-level results")






      