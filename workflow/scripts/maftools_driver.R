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
list.all.maf.files <- lapply(maf.filenames,function(i){
  read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

###merging the all the .maf files
my_maf <- maftools::merge_mafs(list.all.maf.files)
## create out dir
dir.create(file.path(output_path), showWarnings = F)
dir.exists(output_path)
message("Output path is: ",output_path)

## check from maf
aminoacid_cname <- 'Protein_Change'
if('HGVSp' %in% names(my_maf@data)){
  aminoacid_cname <- 'HGVSp'
}
aminoacid_cname

message("Processing MAF for Driver Genes")
dir.create(file.path(output_path,"tables"), showWarnings = F)
dir.create(file.path(output_path,"plots"), showWarnings = F)

tryCatch({
    png(filename = file.path(output_path,"plots","somatic_interactions.png"), width = 500, height = 250, units='mm', res = 400)
    somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.8, colPal = "RdYlGn",
                    sigSymbolsFontSize = 1, sigSymbolsSize = 2.5,showSum = F)
    dev.off()
}, error = function(err) {
    print("Error while plotting Somatic Interaction")
    print(err)
    file.create(file.path(output_path, "plots", "somatic_interactions.png"))
})

tryCatch({
    pdf(file = NULL)
    write.table(as.data.frame(somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.6), row.names = F),
                file = file.path(output_path,"tables","somatic_interactions.tsv"), sep = "\t", row.names = F, col.names = T)
    dev.off()

    my_maf.sig <- oncodrive(maf = my_maf, AACol = aminoacid_cname, minMut = 5, pvalMethod = 'zscore')

    write.table(as.data.frame(my_maf.sig), file = file.path(output_path,"tables","oncodrive.tsv"), sep = "\t", row.names = F, col.names = T)
}, error = function(err) {
    print("Error while creating table")
    print(err)
})

tryCatch({
    png(filename = file.path(output_path,"plots","oncodrive.png"), width = 500, height = 250, units='mm', res = 400)
    plotOncodrive(res = my_maf.sig, fdrCutOff = 0.01, useFraction = T , labelSize = 1)
    dev.off()
}, error = function(err) {
    print("Error while plotting OncoDrive")
    print(err)
    file.create(file.path(output_path, "plots", "oncodrive.png"))
})


tryCatch({
    png(filename = file.path(output_path,"plots","drug_interactions_barplot.png"), width = 500, height = 250, units='mm', res = 400)
    dgi <- drugInteractions(maf = my_maf, fontSize = 0.9, plotType = "bar", top = 25)
    dev.off()
}, error = function(err) {
    print("Error while plotting Drug Interaction")
    print(err)
    file.create(file.path(output_path, "plots", "drug_interactions_barplot.png"))
})


tryCatch({
    png(filename = file.path(output_path,"plots","drug_interactions_piechart.png"), width = 500, height = 250, units='mm', res = 400)
    dgi <- drugInteractions(maf = my_maf, fontSize = 0.7, plotType = "pie", top = 25)
    dev.off()
}, error = function(err) {
    print("Error while plotting Drug Interaction")
    print(err)
    file.create(file.path(output_path, "plots", "drug_interactions_piechart.png"))
})
