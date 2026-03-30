# Create the peak x sample matrix
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

projPath <- '/fh/fast/cohn_l/grp/projects/cuttag-data/CUTTAGpipeline/'

# Read sample sheet to get histone and replicate lists
df <- read.table("data/samples.tsv", header = TRUE, sep = "\t",
                 comment.char = "", strip.white = TRUE)
df$hist <- sub("_.*$", "", df$sample)
df$rep  <- sub("^.*_", "", df$sample)

# save as lists
histL <- unique(df$hist[!df$hist %in% "IgG"])
repL  <- unique(df$rep)

# Master peak list (union across all histones)
mPeak <- GRanges()
# overlap with bam file to get count
for(hist in histL){
  for(rep in repL){
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_control.peaks.stringent.bed"), 
                         header = FALSE, fill = TRUE)
    mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
  }
}
masterPeak = reduce(mPeak)

# BAM-based count matrix
bamDir = paste0(projPath, "alignment/bam")
samples <- paste(rep(histL, each = length(repL)), rep(repL, length(histL)), sep = "_")
countMat = matrix(NA, nrow = length(masterPeak), ncol = length(samples)) # empty matrix of peaks X samples
# set rownames
rownames(countMat) <- paste(seqnames(masterPeak),
                            start(masterPeak),
                            end(masterPeak),
                            sep = ":")
colnames(countMat) <- samples
bamCounts <- countMat

library(chromVAR)
# Overlap with bam file to get count

i = 1
for(hist in histL){
  for(rep in repL){
    bamFile = paste0(bamDir, "/", hist, "_", rep, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}

head(countMat)
# Write raw count matrix (for downstream tools)
write.csv(countMat,
          file.path(projPath, "data/count-matrix/", "raw_count_matrix.csv"),
          row.names = TRUE, quote = FALSE)
cat(nrow(countMat), "peaks in count matrix")

selectR <- rowSums(countMat, na.rm = TRUE) > 5
dataS <- countMat[selectR, , drop = FALSE]
cat(nrow(dataS), "peaks retained")