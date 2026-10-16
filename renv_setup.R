# %%
library(renv)


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# %%
renv::install(c(
"bioc::GenomeInfoDbData",
"bioc::TxDb.Hsapiens.UCSC.hg38.knownGene",
"bioc::BSgenome.Hsapiens.UCSC.hg38",
"bioc::TxDb.Hsapiens.UCSC.hg19.knownGene",
"bioc::BSgenome.Hsapiens.UCSC.hg19"
))


# %%


# Snapshot the environment to lockfile
renv::snapshot()
