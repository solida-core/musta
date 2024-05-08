library("maftools")
library("mclust")
library('NMF')
library('pheatmap')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('stringr')


sessionInfo()
## input
input_maf <- snakemake@input[["mafs"]]
print(input_maf)
project_id <- snakemake@params[["project_id"]]
output_path <- snakemake@params[["outdir"]]
## combine input files
getwd()

maf.filenames <- file.path(input_maf)
list.all.maf.files <- lapply(maf.filenames, function(i) {
    read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

###merging the all the .maf files
my_maf <- maftools::merge_mafs(list.all.maf.files)
## create out dir
dir.create(file.path(output_path), showWarnings = F)
dir.exists(output_path)
message("Output path is: ", output_path)

## check from maf
aminoacid_cname <- 'Protein_Change'
if ('HGVSp' %in% names(my_maf@data)) {
    aminoacid_cname <- 'HGVSp'
}
aminoacid_cname
genome_vers <- 'hg19'
genome_lib <- "BSgenome.Hsapiens.UCSC.hg19"
genome_strings <- c("hg38", "GRCh38")
if (my_maf@data[1, "NCBI_Build"] %in% genome_strings) {
    genome_vers <- 'hg38'
    genome_lib <- "BSgenome.Hsapiens.UCSC.hg38"
    library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
}
genome_vers
genome_lib

chr_prefix <- NULL
add_prefix <- F
# default considera che il maf abbia i chr nominati chr1, chr2 ecc
if (str_detect(my_maf@data[1, "Chromosome"], pattern = "chr", negate = T)) {
    chr_prefix <- 'chr'
    add_prefix <- T
}
chr_prefix
add_prefix

## summaty tables
message("Processing MAF for General Statistics")
dir.create(file.path(output_path, "tables"), showWarnings = F)
dir.create(file.path(output_path, "plots"), showWarnings = F)

## summaty tables
message("Write summary table")
write.table(as.data.frame(my_maf@summary[-c(2), c(1, 2)], row.names = F), file = file.path(output_path, "tables", "overview.tsv"), sep = "\t", row.names = F, col.names = T)

#Shows sample summry.
message("Write Samples Summary table")
write.table(as.data.frame(getSampleSummary(my_maf), row.names = F), file = file.path(output_path, "tables", "sample_summary.tsv"), sep = "\t", row.names = F, col.names = T)

#Shows gene summary.
message("Write Gene Summary table")
write.table(as.data.frame(getGeneSummary(my_maf), row.names = F), file = file.path(output_path, "tables", "gene_summary.tsv"), sep = "\t", row.names = F, col.names = T)

message("Plot MAF Summary")
png(filename = file.path(output_path, "plots", "summary.png"), width = 500, height = 250, units = 'mm', res = 500)
plotmafSummary(maf = my_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T, top = 10, log_scale = F,
               textSize = 4, titleSize = c(1.3, 0.9))
dev.off()


#oncoplot for top ten mutated genes.
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
    'Frame_Shift_Del',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Multi_Hit',
    'Frame_Shift_Ins',
    'In_Frame_Ins',
    'Splice_Site',
    'In_Frame_Del'
)
message("Oncopot generation")
png(filename = file.path(output_path, "plots", "oncoplot.png"), width = 500, height = 250, units = 'mm', res = 400)
oncoplot(maf = my_maf, colors = vc_cols, top = 10, altered = T,
         logColBar = T, drawRowBar = T, draw_titv = F,
         drawBox = F, titleText = NULL, legendFontSize = 1.6, writeMatrix = T, showTumorSampleBarcodes = T, barcode_mar = 15)
dev.off()


png(filename = file.path(output_path, "plots", "large_oncoplot.png"), width = 500, height = 250, units = 'mm', res = 400)
oncoplot(maf = my_maf, colors = vc_cols, altered = F,
         logColBar = F, drawRowBar = F, draw_titv = F,
         drawBox = F, titleText = "Oncoplot 50 genes", legendFontSize = 1.6, drawColBar = F, top = 50, sortByMutation = F)
dev.off()

actual_path <- getwd()
z <- read.table(file.path(actual_path, "onco_matrix.txt"), sep = "\t")

my_maf.titv <- titv(maf = my_maf, plot = FALSE, useSyn = TRUE)
titv(maf = my_maf, plot = FALSE, useSyn = TRUE, file = file.path(output_path, "tables", "transitions_and_transversions"))
#plot titv summary
png(filename = file.path(output_path, "plots", "titv.png"), width = 500, height = 250, units = 'mm', res = 400)
plotTiTv(res = my_maf.titv)
dev.off()


dir.create(file.path(output_path, "plots", "lollipop"), showWarnings = F)
dir.exists(file.path(output_path, "plots", "lollipop"))

for (gene in rownames(z)) {
    print(gene)
    png(filename = file.path(output_path, "plots", "lollipop", paste(gene, ".lollipop.png")), width = 500, height = 250, units = 'mm', res = 400)
    lollipopPlot(
        maf = my_maf,
        gene = gene,
        AACol = aminoacid_cname,
        showMutationRate = TRUE
    )
    dev.off()

}

samples <- as.list(as.data.frame(getSampleSummary(my_maf))["Tumor_Sample_Barcode"])
samples <- as.list(samples$Tumor_Sample_Barcode)

dir.create(file.path(output_path, "plots", "rainfall"), showWarnings = F)
dir.exists(file.path(output_path, "plots", "rainfall"))
i <- 0
for (sample in samples) {
    print(sample)
    i <- i + 1
    print(paste("samples processed:", i, "/", length(samples)))
    try(expr = { png(filename = file.path(output_path, "plots", "rainfall", paste(sample, ".rainfall.png")), width = 500, height = 250, units = 'mm', res = 400)
        rainfallPlot(maf = my_maf, detectChangePoints = TRUE, pointSize = 0.4, ref.build = genome_vers, tsb = sample)
        dev.off() }, silent = TRUE
    )
}


png(filename = file.path(output_path, "plots", "TGCA_compare.png"), width = 500, height = 250, units = 'mm', res = 400)
my_maf.mutload <- tcgaCompare(maf = my_maf, cohortName = project_id, logscale = TRUE, capture_size = 50)
dev.off()

write.table(my_maf.mutload$median_mutation_burden, file = file.path(output_path, "tables", "TGCA_median_mutation_burden.tsv"), sep = "\t", row.names = F, col.names = T)

write.table(my_maf.mutload$mutation_burden_perSample, file = file.path(output_path, "tables", "TGCA_mutation_burden_perSample.tsv"), sep = "\t", row.names = F, col.names = T)

write.table(my_maf.mutload$pairwise_t_test, file = file.path(output_path, "tables", "TGCA_pairwise_t_test.tsv"), sep = "\t", row.names = F, col.names = T)

tryCatch({
    png(filename = file.path(output_path, "plots", "top10_VAF.png"), width = 500, height = 250, units = 'mm', res = 400)
    plotVaf(maf = my_maf, top = 10, showN = T)
    dev.off()
}, error = function(err) {
    print("Error while plotting VAF (top genes: 10)")
    print(err)
    file.create(file.path(output_path, "plots", "top10_VAF.png"))
})

tryCatch({
    png(filename = file.path(output_path, "plots", "top20_VAF.png"), width = 500, height = 250, units = 'mm', res = 400)
    plotVaf(maf = my_maf, top = 20, showN = T)
    dev.off()
}, error = function(err) {
    print("Error while plotting VAF (top genes: 20)")
    print(err)
    file.create(file.path(output_path, "plots", "top20_VAF.png"))
})
