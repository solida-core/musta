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
message("Processing MAF")
dir.create(file.path(output_path,"tables"), showWarnings = F)
dir.create(file.path(output_path,"plots"), showWarnings = F)

dir.create(file.path(output_path,"plots","clusters"), showWarnings = F)
dir.exists(file.path(output_path,"plots","clusters"))

samples=as.list(as.data.frame(getSampleSummary(my_maf))["Tumor_Sample_Barcode"])
samples=as.list(samples$Tumor_Sample_Barcode)
i=0
for(sample in samples){
  print(sample)
  i=i+1
  print(paste("samples processed:",i,"/",length(samples)))
  try(expr = {png(filename = file.path(output_path,"plots","clusters",paste(sample,"_clusters.png",sep = "")), width = 500, height = 250, units='mm', res = 200)
  inferred_het = inferHeterogeneity(maf = my_maf, tsb=sample)
  plotClusters(clusters = inferred_het)
  dev.off()
  write.table(inferred_het$clusterData, file = file.path(output_path,"tables",paste(sample,"_clusterData.tsv",sep = "")), sep = "\t", row.names = F, col.names = T)
  write.table(inferred_het$clusterMeans, file = file.path(output_path,"tables",paste(sample,"_clusterMeans.tsv",sep = "")), sep = "\t", row.names = F, col.names = T)}, silent = TRUE
  )

}

file.create(file.path(output_path,"tables","successful.tsv"))
