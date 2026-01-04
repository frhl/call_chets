
library(ggplot2)
library(data.table)
library(bravastring) # from devtools::install_github("frhl/bravastring")
setwd("~/Projects/06_call_chets/call_chets/examples/set-based/")

# get annotations
anno <- fread("input/variant_annotations.tsv")
anno$position <- NULL

# get single assoc results
d <- fread("output/results.txt.singleAssoc.txt")
d[, expt.p := get_expected_p(p.value)]
d <- merge(d, anno, by.x="MarkerID", by.y="variant_id")

# null variants
calc_inflation(d[d$annotation == "synonymous"]$p.value)

# plot QQ
ggplot(d, aes(x=-log10(expt.p), y=-log10(p.value), color=is_causal)) +
  geom_point() +
  geom_abline()

