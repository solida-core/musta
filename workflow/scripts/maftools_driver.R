library("maftools")
library("mclust")
library('NMF')
library('pheatmap')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('stringr')


sessionInfo()
## input
input_maf <- snakemake@input[["mafs"]]
message("Input files:", input_maf)

project_id <- snakemake@params[["project_id"]]
output_path <- snakemake@params[["outdir"]]
all_variants <- as.logical(snakemake@params[["all_variants"]])
## combine input files
getwd()

maf.filenames <- file.path(input_maf)
list.all.maf.files <- lapply(maf.filenames,function(i){
  read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

###merging the all the .maf files
maf_all <- maftools::merge_mafs(list.all.maf.files)

if (all_variants == TRUE) {
    my_maf <- maf_all
} else {
    if ("FILTER" %in% getFields(maf_all)) {
        my_maf <- subsetMaf(maf_all, query = "FILTER == 'PASS'")
    } else if ("SOMATIC" %in% getFields(maf_all)) {
        my_maf <- subsetMaf(maf_all, query = "SOMATIC == 'true'")
    } else {
        my_maf <- maf_all
    }
}


## create out dir
dir.create(file.path(output_path), showWarnings = F)
dir.exists(output_path)
message("Output path is: ",output_path)

## check from maf
aminoacid_cname <- 'Protein_Change'
if('HGVSp' %in% names(my_maf@data)){
  aminoacid_cname <- 'HGVSp'
}
# aminoacid_cname

message("Processing MAF for Driver Genes")
dir.create(file.path(output_path,"tables"), showWarnings = F)
dir.create(file.path(output_path,"plots"), showWarnings = F)


message("Somatic Interation: detect mutually exclusive, co-occuring and altered genesets.")
tryCatch({
    png(filename = file.path(output_path,"plots","somatic_interactions.png"), width = 500, height = 250, units='mm', res = 400)
    somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.8, colPal = "RdYlGn",
                    sigSymbolsFontSize = 1, sigSymbolsSize = 2.5,showSum = F)
    dev.off()

    pdf(file = NULL)
    write.table(as.data.frame(somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.6), row.names = F),
                file = file.path(output_path,"tables","somatic_interactions.tsv"), sep = "\t", row.names = F, col.names = T)
    dev.off()

}, error = function(err) {
    print("Error while processing Somatic Interaction")
    print(err)
})


message("Checks for drug-gene interactions and druggable categories ")
tryCatch({
    png(filename = file.path(output_path,"plots","drug_interactions.barplot.png"), width = 500, height = 250, units='mm', res = 400)
    dgi <- drugInteractions(maf = my_maf, fontSize = 0.9, plotType = "bar", top = 25)
    dev.off()

    png(filename = file.path(output_path,"plots","drug_interactions.piechart.png"), width = 500, height = 250, units='mm', res = 400)
    dgi <- drugInteractions(maf = my_maf, fontSize = 0.7, plotType = "pie", top = 25)
    dev.off()
}, error = function(err) {
    print("Error while checking for drug-gene interactions and druggable categories ")
    print(err)
})


message("OncoDrive: Detect cancer driver genes based on positional clustering of variants.")
tryCatch({
    my_maf.sig <- oncodrive(maf = my_maf, AACol = aminoacid_cname, minMut = 5, pvalMethod = 'zscore')
    write.table(as.data.frame(my_maf.sig), file = file.path(output_path,"tables","oncodrive.tsv"), sep = "\t", row.names = F, col.names = T)

    png(filename = file.path(output_path,"plots","oncodrive.png"), width = 500, height = 250, units='mm', res = 400)
    plotOncodrive(res = my_maf.sig, fdrCutOff = 0.01, useFraction = T , labelSize = 1)
    dev.off()
}, error = function(err) {
    print("Error while detecting cancer driver genes")
    print(err)
})
