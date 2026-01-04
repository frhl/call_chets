
library(ggplot2)
library(data.table)
library(bravastring) # from devtools::install_github("frhl/bravastring")
setwd("~/Projects/06_call_chets/call_chets/examples/set-based/")

# get annotations
anno <- fread("input/variant_annotations.tsv")
anno$position <- NULL

# get variant-level results
d <- fread("output/saige.step2.additive.variant.txt")
d[, expt.p := get_expected_p(p.value)]
d <- merge(d, anno, by.x="MarkerID", by.y="variant_id")

# plot QQ
ggplot(d, aes(x=-log10(expt.p), y=-log10(p.value), color=is_causal)) +
  geom_point() +
  geom_abline()

# get gene-level results
d <- fread("output/saige.step2.additive.group.txt")
d[, expt.p := get_expected_p(Pvalue)]

ggplot(d, aes(x=-log10(expt.p), y=-log10(Pvalue), color=Group)) +
  geom_point() +
  geom_abline()




