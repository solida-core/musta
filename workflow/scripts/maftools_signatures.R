library("maftools")
library("mclust")
library('NMF')
library('pheatmap')
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library('stringr')
library("tools")

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
list.all.maf.files <- lapply(maf.filenames, function(i) {
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
message("Output path is: ", output_path)

## check from maf
aminoacid_cname <- 'Protein_Change'
if ('HGVSp' %in% names(my_maf@data)) {
    aminoacid_cname <- 'HGVSp'
}
# aminoacid_cname

genome_vers <- 'hg19'
genome_lib <- "BSgenome.Hsapiens.UCSC.hg19"
genome_strings <- c("hg38", "GRCh38")
if (my_maf@data[1, "NCBI_Build"] %in% genome_strings) {
    genome_vers <- 'hg38'
    genome_lib <- "BSgenome.Hsapiens.UCSC.hg38"
    library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
}
message("Genome Version: ", genome_vers)
message("Genome Library: ", genome_lib)

chr_prefix <- NULL
add_prefix <- F
# default considera che il maf abbia i chr nominati chr1, chr2 ecc
if (str_detect(my_maf@data[1, "Chromosome"], pattern = "chr", negate = T)) {
    chr_prefix <- 'chr'
    add_prefix <- T
}

# chr_prefix
# add_prefix

## summaty tables
message("Processing MAF for Signatures")
dir.create(file.path(output_path, "tables"), showWarnings = F)
dir.create(file.path(output_path, "plots"), showWarnings = F)

message("Extract single 5' and 3' bases flanking the mutated site for de-novo signature analysis.")
tryCatch({
    my_maf.tnm <- trinucleotideMatrix(maf = my_maf, prefix = chr_prefix, add = add_prefix, ref_genome = genome_lib)

    nmf_matrix <- as.data.frame(my_maf.tnm$nmf_matrix)
    nmf_matrix <- cbind(rownames(nmf_matrix), nmf_matrix)
    names(nmf_matrix)[names(nmf_matrix) == "rownames(nmf_matrix)"] <- "sample_ID"
    write.table(nmf_matrix, file = file.path(output_path, "tables", "nmf_matrix.tsv"), sep = "\t", row.names = F, col.names = T)
}, error = function(err) {
    print("Error while extracting single 5' and 3' bases flanking the mutated site for de-novo signature analysis.")
    print(err)
})

message("Estimates APOBEC enrichment scores.")
tryCatch({
    APOBEC_scores_matrix <- as.data.frame(my_maf.tnm$APOBEC_scores)
    write.table(APOBEC_scores_matrix, file = file.path(output_path, "tables", "APOBEC_scores_matrix.tsv"), sep = "\t", row.names = F, col.names = T)
}, error = function(err) {
    print("Error while Estimating APOBEC enrichment scores")
    print(err)
})

message("Plot differences between APOBEC enriched and non-APOBEC enriched samples.")
tryCatch({
    apobec_diff_df <- plotApobecDiff(tnm = my_maf.tnm, maf = my_maf, pVal = 0.05)
    apobec_diff_df$SampleSummary

    png(filename = file.path(output_path, "plots", "Apobecdiff.png"), width = 500, height = 250, units = 'mm', res = 400)
    plotApobecDiff(tnm = my_maf.tnm, maf = my_maf, pVal = 0.2)
    dev.off()
}, error = function(err) {
    print("Error while plotting differences between APOBEC enriched and non-APOBEC enriched samples.")
    print(err)
})


message("Estimate number of signatures")
tryCatch({
    my_maf.sign <- estimateSignatures(mat = my_maf.tnm, nTry = 6)
}, error = function(err) {
    my_maf.sign <- estimateSignatures(mat = my_maf.tnm, nTry = 6,  pConstant=0.1)
    # print("Error while extracting Signatures")
    # print(err)
})

message("Draw an elbow plot of cophenetic correlation metric.")
tryCatch({
    png(filename = file.path(output_path, "plots", "plotCophenetic.png"), width = 500, height = 250, units = 'mm', res = 400)
    plotCophenetic(res = my_maf.sign)
    dev.off()
}, error = function(err) {
    print("Error while drawing an elbow plot of cophenetic correlation metric. ")
    print(err)
    # file.create(file.path(output_path, "plots", "plotCophenetic.png"))
})

tryCatch({
    new_df <- as.data.frame(my_maf.sign$nmfSummary$cophenetic)
    names(new_df) <- "cophenetic"
    rownames(new_df) <- rownames(my_maf.sign$nmfSummary)
    new_df["rank"] <- rownames(my_maf.sign$nmfSummary)
    new_df["drop"] <- NA
    new_df[1, "drop"] <- 0
    for (row in 2:nrow(new_df)) {
        new_df[row, "drop"] <- new_df[row, "cophenetic"] - new_df[row - 1, "cophenetic"]
    }
}, error = function(err) {
    print("Error")
    print(err)
})

message("Extract mutational signatures")
tryCatch({
    signature_to_extract <- as.integer(new_df[new_df$drop == min(new_df$drop), "rank"])

    tryCatch({
        my_maf.sig <- extractSignatures(mat = my_maf.tnm, n = signature_to_extract)
    }, error = function(err) {
        my_maf.sig <- extractSignatures(mat = my_maf.tnm, n = signature_to_extract, pConstant=0.1)
    }

    signatures_matrix <- as.data.frame(my_maf.sig$signatures)
    signatures_matrix <- cbind(rownames(signatures_matrix), signatures_matrix)
    names(signatures_matrix)[names(signatures_matrix) == "rownames(signatures_matrix)"] <- "trinucleotide"
    write.table(signatures_matrix, file = file.path(output_path, "tables", "signatures_matrix.tsv"), sep = "\t", row.names = F, col.names = T)

    message("Compate against updated version3 60 signatures")
    my_maf.v3.cosm <- compareSignatures(nmfRes = my_maf.sig, sig_db = "SBS")

    cosine_simil_matrix <- as.data.frame(my_maf.v3.cosm$cosine_similarities)
    cosine_simil_matrix <- cbind(rownames(cosine_simil_matrix), cosine_simil_matrix)
    names(cosine_simil_matrix)[names(cosine_simil_matrix) == "rownames(cosine_simil_matrix)"] <- "signature"
    write.table(cosine_simil_matrix, file = file.path(output_path, "tables", "cosime_similarities_matrix.tsv"), sep = "\t", row.names = F, col.names = T)

    dev.off()

    png(filename = file.path(output_path, "plots", "cosine_similarities_heatmap.png"), width = 300, height = 150, units = 'mm', res = 400)
    pheatmap(mat = my_maf.v3.cosm$cosine_similarities, cluster_rows = FALSE, main = "Cosine similarity against validated signatures",
             fontsize_col = 8, cellwidth = 10, cellheight = 10, fontsize_row = 8)
    dev.off()

    png(filename = file.path(output_path, "plots", "cosmic_signatures.png"), width = 500, height = 250, units = 'mm', res = 400)
    plotSignatures(nmfRes = my_maf.sig, title_size = 1.5, sig_db = "SBS", font_size = 1.4, yaxisLim = NA)
    dev.off()

    png(filename = file.path(output_path, "plots", "signature_contributions.png"), width = 500, height = 250, units = 'mm', res = 400)
    plotSignatures(nmfRes = my_maf.sig, title_size = 1.2, sig_db = "SBS", contributions = T, show_barcodes = F)
    dev.off()
}, error = function(err) {
    print("Error while extracting Mutational Signatures")
    print(err)
    file.create(file.path(output_path, "plots", "cosmic_signatures.png"))
})
