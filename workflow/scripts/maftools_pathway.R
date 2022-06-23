library("maftools")
library("mclust")
library('NMF')
library('pheatmap')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('stringr')


sessionInfo()
## input
input_maf=snakemake@input[["mafs"]]
print(input_maf)
project_id=snakemake@params[["project_id"]]
output_path=snakemake@params[["outdir"]]
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
aminoacid_cname='Protein_Change'
if('HGVSp' %in% names(my_maf@data)){
  aminoacid_cname='HGVSp'
}
aminoacid_cname
genome_vers='hg19'
genome_lib="BSgenome.Hsapiens.UCSC.hg19"
genome_strings=c("hg38","GRCh38")
if(my_maf@data[1,"NCBI_Build"] %in% genome_strings){
  genome_vers='hg38'
  genome_lib="BSgenome.Hsapiens.UCSC.hg38"
  library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
}
genome_vers
genome_lib

chr_prefix=NULL
add_prefix=F
# default considera che il maf abbia i chr nominati chr1, chr2 ecc
if(str_detect(my_maf@data[1,"Chromosome"], pattern = "chr", negate = T)){
  chr_prefix='chr'
  add_prefix=T
}
chr_prefix
add_prefix

## summaty tables
message("Processing MAF for General Statistics")
dir.create(file.path(output_path,"tables"), showWarnings = F)
dir.create(file.path(output_path,"plots"), showWarnings = F)

png(filename = file.path(output_path,"plots","oncogenic_pathways.png"), width = 500, height = 250, units='mm', res = 400)
OncogenicPathways(maf = my_maf, fontSize = 1.5)
dev.off()

pathways=as.data.frame(OncogenicPathways(maf = my_maf))
write.table(pathways, file = file.path(output_path,"tables","oncogenic_pathways.tsv"), sep = "\t", row.names = F, col.names = T)

dir.create(file.path(output_path,"plots","pathways"), showWarnings = F)
dir.exists(file.path(output_path,"plots","pathways"))
for(pathway in pathways$Pathway){
  png(filename = file.path(output_path,"plots","pathways",paste(pathway,"_pathway.png",sep = "")), width = 500, height = 250, units='mm', res = 400)
  PlotOncogenicPathways(maf = my_maf, pathways = pathway, showTumorSampleBarcodes = F, SampleNamefontSize = 0.4)
  dev.off()
}
