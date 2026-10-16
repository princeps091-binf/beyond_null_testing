# %%

library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)

# %%

ah <- AnnotationHub()

# %%

# Define promoters as TSS +/- 1000 bp for standard human protein-coding genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoters_gr <- promoters(genes(txdb), upstream = 1000, downstream = 1000)
promoters_gr <- keepStandardChromosomes(promoters_gr, pruning.mode = "coarse")

# %%

# Keep a robust subset (e.g., chromosome 1 & 2 for rapid download and extraction)
promoters_gr <- promoters_gr[seqnames(promoters_gr) %in% c(paste0('chr',1:22),'chrX','chrY')]

message("--> 2. Calculating exact GC content from BSgenome...")
promoter_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, promoters_gr)
promoters_gr$gc_content <- letterFrequency(promoter_seqs, letters = "GC", as.prob = TRUE)[, 1]

# %%
# Query real ENCODE K562 POLR2A ChIP-seq and control tracks
active_mark <- query(ah, c("K562",'E123',"H3K4me3",'fc'))
dnase_peak <- query(ah, c("K562","E123","dnase",'narrow'))

# Extract resources (taking the top matching records)
active_mark_bw   <- ah[[names(active_mark)]]

dnase_peak   <- ah[[names(dnase_peak)]]

# %%
active_scores <- import(active_mark_bw, which = promoters_gr, format = "BigWig")
active_cov <- coverage(active_scores, weight = "score")
# %%
promoters_gr$active_ip_signal <- pmap_dbl(
  list(
    chr = as.character(seqnames(promoters_gr)),
    start = start(promoters_gr),
    end = end(promoters_gr)
  ),
  function(chr, start, end) {
    if (!chr %in% names(active_cov)) return(0)
    
    chr_len <- length(active_cov[[chr]])
    
    # 2. Clamp coordinates to valid chromosome boundaries
    start_safe <- max(1, start)
    end_safe <- min(chr_len, end)
    
    # 3. Handle edge case where the range is entirely out of bounds
    if (start_safe > end_safe) return(0)
    
    # 4. Extract window safely and compute max
    window_signal <- active_cov[[chr]][start_safe:end_safe]
    if (length(window_signal) == 0 || all(is.na(window_signal))) return(0)
    
    max(window_signal, na.rm = TRUE)
  }
)
# %%
# C. ATAC-seq peak overlap (binary indicator)
hits <- findOverlaps(promoters_gr, dnase_peak)
promoters_gr$dnase_peak <- 0
promoters_gr$dnase_peak[unique(queryHits(hits))] <- 1

# %%

encode_df <- as.data.frame(promoters_gr) %>%
  as_tibble() %>%
  rename(promoter_id = gene_id) %>%
  filter(!is.na(active_ip_signal), !is.na(gc_content)) %>%
  select(promoter_id, active_ip_signal, dnase_peak, gc_content) 

write_csv(encode_df, "./data/encode_k562_promoters.csv")


