library("maftools")
library("mclust")
library('NMF')
library('pheatmap')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('stringr')


sessionInfo()
## input
input_maf=snakemake@input[["mafs"]]

#project_id="project_id"
#output_path="/home/matteo/Scrivania/aleale_new_20"
#input_maf="test-data-somatic/data/maf/BRCA_subset_25_samples.maf"

print(input_maf)
output_path=snakemake@params[["outdir"]]
## combine input files
dirname(input_maf)
setwd(dirname(input_maf[1]))
getwd()
###reading mutiple .maf files as a large list

#maf.filenames <- list.files(full.names=TRUE, pattern = ".maf")
maf.filenames <- file.path(".",basename(input_maf))
list.all.maf.files <- lapply(maf.filenames,function(i){
  read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

###merging the all the .maf files
my_maf <- maftools::merge_mafs(list.all.maf.files)
#write.table(my_maf, file = file.path(output_path,"merged_maf.tsv"), sep = "\t", row.names = F, col.names = T)
## create out dir
dir.create(file.path(output_path), showWarnings = F)
dir.exists(output_path)
setwd(output_path)
message("Output path is: ",output_path)

## read maf
#my_maf = read.maf(maf = input_maf)

## check from maf
aminoacid_cname='Protein_Change'
if('HGVSp' %in% names(my_maf@data)){
  aminoacid_cname='HGVSp'
}
aminoacid_cname

message("Processing MAF for Driver Genes")
dir.create(file.path(output_path,"tables"), showWarnings = F)
dir.create(file.path(output_path,"plots"), showWarnings = F)

png(filename = file.path(output_path,"plots","somatic_interactions.png"), width = 500, height = 250, units='mm', res = 400)
somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.8, colPal = "RdYlGn",
                    sigSymbolsFontSize = 1, sigSymbolsSize = 2.5,showSum = F)
dev.off()

pdf(file = NULL)
write.table(as.data.frame(somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.6), row.names = F),
            file = file.path(output_path,"tables","somatic_interactions.tsv"), sep = "\t", row.names = F, col.names = T)
dev.off()

my_maf.sig = oncodrive(maf = my_maf, AACol = aminoacid_cname, minMut = 5, pvalMethod = 'zscore')

write.table(as.data.frame(my_maf.sig), file = file.path(output_path,"tables","oncodrive.tsv"), sep = "\t", row.names = F, col.names = T)

png(filename = file.path(output_path,"plots","oncodrive.png"), width = 500, height = 250, units='mm', res = 400)
plotOncodrive(res = my_maf.sig, fdrCutOff = 0.01, useFraction = T , labelSize = 1)
dev.off()

png(filename = file.path(output_path,"plots","drug_interactions_barplot.png"), width = 500, height = 250, units='mm', res = 400)
dgi = drugInteractions(maf = my_maf, fontSize = 0.9,plotType = "bar",top = 25)
dev.off()

png(filename = file.path(output_path,"plots","drug_interactions_piechart.png"), width = 500, height = 250, units='mm', res = 400)
dgi = drugInteractions(maf = my_maf, fontSize = 0.7,plotType = "pie", top = 25)
dev.off()
