library("maftools")
library("mclust")
library('NMF')
library('pheatmap')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('stringr')


## input
input_maf=snakemake@input[[1]]
project_id=snakemake@params[["project_id"]]
output_path=snakemake@params[["outdir"]]
## combine input files
dirname(input_maf)
setwd(dirname(input_maf[1]))
###reading mutiple .maf files as a large list
maf.filenames <- list.files(full.names=TRUE, pattern = "_funcotated.maf")
list.all.maf.files <- lapply(maf.filenames,function(i){
  read.delim(i, sep = "\t", header = TRUE, fill = TRUE, comment.char = "#")
})

###merging the all the .maf files
my_maf <- maftools::merge_mafs(list.all.maf.files)

## create out dir
dir.create(file.path(output_path), showWarnings = F)
dir.exists(output_path)
setwd(output_path)

## read maf
#my_maf = read.maf(maf = input_maf)

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
write.table(as.data.frame(my_maf@summary[-c(2),c(1,2)], row.names = F), file = file.path(output_path,"overview.tsv"), sep = "\t", row.names = F, col.names = T)

list(my_maf@clinical.data)

#Shows sample summry.
write.table(as.data.frame(getSampleSummary(niasmic), row.names = F), file = file.path(output_path,"sample_summary.tsv"), sep = "\t", row.names = F, col.names = T)

#Shows gene summary.
write.table(as.data.frame(getGeneSummary(niasmic), row.names = F), file = file.path(output_path,"gene_summary.tsv"), sep = "\t", row.names = F, col.names = T)

plotmafSummary(maf = my_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T, top = 10, log_scale = F)
dev.off()
png(filename = file.path(output_path,"summary.png"), width = 500, height = 250, units='mm', res = 500)
plotmafSummary(maf = my_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T, top = 10, log_scale = T)
dev.off()


#oncoplot for top ten mutated genes.
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

png(filename = file.path(output_path,"oncoplot.png"), width = 500, height = 250, units='mm', res = 400)
oncoplot(maf = my_maf, colors = vc_cols, top = 10, altered = T, 
         logColBar = T, drawRowBar = T,draw_titv = F, 
         drawBox = F, titleText = NULL, legendFontSize = 1.6, writeMatrix = T, showTumorSampleBarcodes = T, barcode_mar = 15)
dev.off()


png(filename = file.path(output_path,"large_oncoplot.png"), width = 500, height = 250, units='mm', res = 400)
oncoplot(maf = my_maf, colors = vc_cols, altered = F, 
         logColBar = F, drawRowBar = F,draw_titv = F, 
         drawBox = F, titleText = "Oncoplot 50 genes", legendFontSize = 1.6, drawColBar = F, top = 50, sortByMutation = F)
dev.off()


z=read.table(file.path(output_path,"onco_matrix.txt"), sep = "\t")
rownames(z)


my_maf.titv = titv(maf = my_maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
png(filename = file.path(output_path,"titv.png"), width = 500, height = 250, units='mm', res = 400)
plotTiTv(res = my_maf.titv)
dev.off()

dir.create(file.path("lollipop"), showWarnings = F)
dir.exists("lollipop")

for(gene in rownames(z)){
  print(gene)
  png(filename = file.path(output_path,"lollipop",paste(gene,".lollipop.png")), width = 500, height = 250, units='mm', res = 400)
  lollipopPlot(
    maf = my_maf,
    gene = gene,
    AACol = aminoacid_cname,
    showMutationRate = TRUE,
  )
  dev.off()
  
}

samples=as.list(my_maf@clinical.data$Tumor_Sample_Barcode)

dir.create(file.path("rainfall"), showWarnings = F)
dir.exists("rainfall")

for(sample in samples){
  print(sample)
  png(filename = file.path(output_path,"rainfall",paste(sample,".rainfall.png")), width = 500, height = 250, units='mm', res = 400)
  rainfallPlot(maf = my_maf, detectChangePoints = TRUE, pointSize = 0.4, ref.build = genome_vers, tsb=sample)
  dev.off()
}


png(filename = file.path(output_path,"TGCA_compare.png"), width = 500, height = 250, units='mm', res = 400)
my_maf.mutload = tcgaCompare(maf = my_maf, cohortName = project_id, logscale = TRUE, capture_size = 50)
dev.off()


png(filename = file.path(output_path,"top10_VAF.png"), width = 500, height = 250, units='mm', res = 400)
plotVaf(maf = my_maf, top = 10, showN = T)
dev.off()

png(filename = file.path(output_path,"top20_VAF.png"), width = 500, height = 250, units='mm', res = 400)
plotVaf(maf = my_maf, top = 20, showN = T)
dev.off()



png(filename = file.path(output_path,"somatic_interactions.png"), width = 500, height = 250, units='mm', res = 400)
somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.6)
dev.off()

write.table(as.data.frame(somaticInteractions(maf = my_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.6), row.names = F), 
            file = file.path(output_path,"somatic_interactions.tsv"), sep = "\t", row.names = F, col.names = T)



my_maf.sig = oncodrive(maf = my_maf, AACol = aminoacid_cname, minMut = 5, pvalMethod = 'zscore')
#niasmic.sig = oncodrive(maf = niasmic, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
head(my_maf.sig)

write.table(my_maf.sig, file = file.path(output_path,"oncodrive.tsv"), sep = "\t", row.names = F, col.names = T)

png(filename = file.path(output_path,"oncodrive.png"), width = 500, height = 250, units='mm', res = 400)
plotOncodrive(res = my_maf.sig, fdrCutOff = 0.01, useFraction = TRUE, labelSize = 0.8)
dev.off()

png(filename = file.path(output_path,"drug_interactions.png"), width = 500, height = 250, units='mm', res = 400)
dgi = drugInteractions(maf = my_maf, fontSize = 0.9)
dev.off()


png(filename = file.path(output_path,"oncogenic_pathways.png"), width = 500, height = 250, units='mm', res = 400)
OncogenicPathways(maf = my_maf)
dev.off()

nsa=OncogenicPathways(maf = my_maf)
write.table(nsa, file = file.path(output_path,"oncogenic_pathways.tsv"), sep = "\t", row.names = F, col.names = T)

dir.create(file.path("pathways"), showWarnings = F)
dir.exists("pathways")
for(pathway in nsa$Pathway){
  numsamples=nsa
  png(filename = file.path(output_path,"pathways",paste(pathway,"_pathway.png")), width = 500, height = 250, units='mm', res = 400)
  PlotOncogenicPathways(maf = my_maf, pathways = pathway, showTumorSampleBarcodes = F, SampleNamefontSize = 0.4)
  dev.off()
  ## sample font size può essere variabile in base al numero di cmapioni
}


dir.create(file.path("clusters"), showWarnings = F)
dir.exists("clusters")
for(sample in samples){
  print(sample)
  png(filename = file.path(output_path,"pathways",paste(sample,"_clusters.png")), width = 500, height = 250, units='mm', res = 200)
  i = inferHeterogeneity(maf = my_maf, tsb=sample)
  plotClusters(clusters = i)
  dev.off()
}




my_maf.tnm = trinucleotideMatrix(maf = my_maf, prefix = chr_prefix, add = add_prefix, ref_genome = genome_lib)


png(filename = file.path(output_path,"Apobecdiff.png"), width = 500, height = 250, units='mm', res = 400)
plotApobecDiff(tnm = my_maf.tnm, maf = my_maf, pVal = 0.2)
dev.off()


my_maf.sign = estimateSignatures(mat = my_maf.tnm, nTry = 6)

plotCophenetic(res = my_maf.sign)

my_maf.sig = extractSignatures(mat = my_maf.tnm, n = 3)

#Compate against updated version3 60 signatures 
my_maf.v3.cosm = compareSignatures(nmfRes = my_maf.sig, sig_db = "SBS")

pheatmap::pheatmap(mat = my_maf.v3.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")


png(filename = file.path(output_path,"cosmic_signatures.png"), width = 500, height = 250, units='mm', res = 400)
plotSignatures(nmfRes = my_maf.sig, title_size = 1.2, sig_db = "SBS")
dev.off()

png(filename = file.path(output_path,"signature_contributions.png"), width = 500, height = 250, units='mm', res = 400)
plotSignatures(nmfRes = my_maf.sig, title_size = 1.2, sig_db = "SBS", contributions = T, show_barcodes = T)
dev.off()


#### end
