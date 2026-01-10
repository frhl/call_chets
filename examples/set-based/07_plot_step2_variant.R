# plot vairant level results

library(ggplot2)
library(data.table)
library(bravastring) # from devtools::install_github("frhl/bravastring")
setwd("~/Projects/06_call_chets/call_chets/examples/set-based/")

# get non-additive (dominance deviation) results
cols_to_keep <- c("MarkerID", "BETA", "SE", "Tstat", "var", "p.value")
d_nonadd <- fread("output/saige.step2.nonadditive.variant.txt")
d_nonadd <- d_nonadd[,..cols_to_keep]

# merge with additive results
d_add <- fread("output/saige.step2.additive.variant.txt")
d <- merge(d_add, d_nonadd, by="MarkerID", suffixes = c(".add", ".dom"))

# let's also add recessive results
d_rec <- fread("output/saige.step2.recessive.variant.txt")
d_rec <- d_rec[,..cols_to_keep]
colnames(d_rec)[2:6] <- paste0(colnames(d_rec)[2:6], ".rec")
d <- merge(d, d_rec, by="MarkerID")

# merge with annotations
anno <- fread("input/variant_annotations.tsv")
d <- merge(d, anno, by.x="MarkerID", by.y="variant_id")

# restrict variants with enough homozygotes alternates
d <- d[(n_hom_alt>=5 & n_hom_alt<50) & maf < 0.05]

# calculate expected P
d[, p.value.expt.add := get_expected_p(p.value.add)]
d[, p.value.expt.dom := get_expected_p(p.value.dom)]
d[, p.value.expt.rec := get_expected_p(p.value.rec)]

# We also add a 'model' column to distinguish between them
plot_data <- rbind(
  d[, .(p.value = p.value.add, expt.p = p.value.expt.add, is_causal, model = "Additive [0,1,2]")],
  d[, .(p.value = p.value.dom, expt.p = p.value.expt.dom, is_causal, model = "Non-additive")],
  d[, .(p.value = p.value.rec, expt.p = p.value.expt.rec, is_causal, model = "Recessive [0,0,2]")]
)

# 2. Plot side-by-side
ggplot(plot_data, aes(x = -log10(expt.p), y = -log10(p.value), color = is_causal)) +
  geom_point(data=plot_data[is_causal==FALSE], size=2.5) +
  geom_point(data=plot_data[is_causal==TRUE], size=2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("FALSE" = "#999999", "TRUE" = "#E69F00")) +
  facet_wrap(~model) + 
  theme_bw() +
  labs(
    x = "-log10(Expected P-value)",
    y = "-log10(Observed P-value)",
    color = "Is Causal"
  ) + 
  ggtitle("Simulating a completely recessive trait", "variant-level results")




