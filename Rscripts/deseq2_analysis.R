###########################################################################
###########################################################################
###
### DIFFERENTIAL GENE EXPRESSION ANALYSIS USING DESEQ2
###
###########################################################################
###########################################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library(here)
source(here("Rscripts/load_libraries.R"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
install.packages("IHW")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ANALYSIS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................DESeq2 pipeline/including featurecounts

#................................genome annotation file
gtf_file <- here("data/genome_data/CP023154.gtf")

#................................mapped files
file_folder <- list.files(path = here("data/mapped_data/star"), full.names = T)
file_path <- list.files(path = here("data/mapped_data/star"), full.names = F)
bam_list <- str_c(file_folder,"/", file_path, "_Aligned.sortedByCoord.out.bam",sep = "") 

#................................calculate count matrix using featurecounts and label features
counts <- featureCounts(bam_list[1:6],verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "exon", allowMultiOverlap = F, isLongRead = F)$counts
colnames_data <- c("replicate", "strain", "treatment")
colData <- matrix(ncol = 3, nrow = 6)
colnames(colData) <- colnames_data 
colData <- data.table(colData) %>%
  mutate(replicate = as.factor(rep(1:3,2)),
         strain = as.factor(c(rep(52,6))),
         treatment = as.factor(c(rep("no",3), rep("yes",3))))
counts_table <- as_tibble(counts)
counts_table$gene <- rownames(counts)

#................................DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                      colData = colData,
                                      design = ~treatment)

#................................remove genes with no counts
dds <- dds[ rowSums(counts(dds)) > 1, ]

#................................transformation of data / vst (variance stabilizing transformation for negative binomial data with a dispersion-mean trend)
vsd <- vst(dds, blind = T)
results_table <- as_tibble(assay(vsd)) 
colnames(results_table) <- as.factor(1:6) 
table <- results_table %>%  
  mutate(means_ = rowMeans(.[1:6]),
         sd_ = rowSds(as.matrix(.[1:6])),
         genename = rownames(.))

#................................PCA PLOT (FIG. 3a)
pcaData <- plotPCA(vsd, intgroup = c("replicate", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(here("/figures/pca_plot.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment, shape = replicate, fill = treatment)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_shape_manual(values = c(22,23,24)) +
  coord_fixed(ratio = 3) +
  theme_light() +
  theme_Publication_white() +
  scale_color_viridis(discrete = TRUE, option = "magma", direction = 1, end = 0.75, begin = 0.1) +
  scale_fill_viridis(discrete = TRUE, option = "magma", direction = 1, end = 0.75, begin = 0.1)
dev.off()

#................................BEESWARM PLOT FOR PF0740 (FIG. 3b)
coppertransporter <- rownames(vsd)[which(rownames(vsd) == "PFDSM3638_03695")]
geneCounts <- plotCounts(dds, gene = coppertransporter, intgroup = c("replicate","treatment"), returnData = T)

pdf(here("figures/CopA_counts.pdf"), 
    width = 5, height = 7, paper = "special",onefile=FALSE)
ggplot(geneCounts, aes(x = treatment, y = as.numeric(count), color = treatment, shape = replicate, fill = treatment)) +
  scale_y_log10(limits = c(90,10000)) +  
  geom_beeswarm(cex = 5, size = 5, alpha = 1, dodge.width = 0.3) +
  scale_shape_manual(values = c(22,23,24)) +
  ylab("counts") +
  xlab("PF0740: copper transporter") +
  theme_Publication_white() +
  scale_color_viridis(discrete = TRUE, option = "magma", direction = 1, end = 0.75, begin = 0.1) +
  scale_fill_viridis(discrete = TRUE, option = "magma", direction = 1, end = 0.75, begin = 0.1)
dev.off()

#................................calculate differential gene expression and add significance information to genes / apeglm transformation (apeglm provides empirical Bayes shrinkage estimators for effect sizes for a variety of GLM models; apeglm stands for “Approximate Posterior Estimation for GLM”. apeglm package version: 1.7.5)
dds_results <- DESeq(dds)
res <- lfcShrink(dds_results, coef="treatment_yes_vs_no", type="apeglm")

#................................datapoint with NA of padj have too less counts and are therefore excluded
res_ggplot <- res %>%
  as_tibble() %>%
  rownames_to_column() %>%
  mutate(padj = ifelse(is.na(padj) == T, 1,padj)) %>%
  mutate(type_of_regulation = ifelse(log2FoldChange > 1 & padj < 0.05, "up", 
                                     ifelse(log2FoldChange < -1 & padj < 0.05, "down",
                                            ifelse(is.na(padj), "up","not significant"))),
         id = rownames(res))



#................................MA PLOT OF ALL GENES THAT ARE REGULATED UPON COPPER SHOCK
pdf(here("/figures/ma_plot.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggplot(data = res_ggplot, aes(x = baseMean, y = log2FoldChange, color = type_of_regulation, fill = type_of_regulation,text = rowname)) +
  geom_point(size = 4, stroke = 1.5, shape = 21) +
  scale_color_viridis(option = "magma", discrete = T, begin = 0.1, end = 0.6, alpha = 1) +
  scale_fill_viridis(option = "magma", discrete = T, begin = 0.1, end = 0.6, alpha = 0.6) +
  scale_x_log10(limits = c(0.1,1e6)) +
  scale_y_continuous(limits = c(-8,8)) +
  xlab("mean counts") +
  ylab("log2FoldChange") +
  theme_Publication_white() +
  geom_hline(yintercept = 0, color = viridis_pal(option = "magma")(10)[1], alpha = 0.5, size = 1) +
  geom_hline(yintercept = 1, color = viridis_pal(option = "magma")(10)[1], alpha = 0.5, size = 1, linetype = 2) +
  geom_hline(yintercept = -1, color = viridis_pal(option = "magma")(10)[1], alpha = 0.5, size = 1, linetype = 2)
dev.off()

#................................add arcog information and save to file
arcog_pfu <- fread(here("data/arcog_data/arcog_pfu.tsv"))
res_ggplot_arcog <- res_ggplot %>%
  left_join(arcog_pfu, by = c("id" = "gene"))
write.table(res_ggplot_arcog, file = here("tables/deseq2_table.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")


#................................HEATMAP ANALYSIS
#..............................find 25 most-significantly regulated genes
top25  <- rownames(res)[order(res$log2FoldChange, decreasing = T)[1:25]]

# > alternative plot: Intensity (log10 counts values) AND Fold changes
count_heatmap <- counts(dds_results) %>%
  as_tibble()
colnames(count_heatmap) <- c("52_1", "52_2", "52_3", "52_1_cu", "52_2_cu", "52_3_cu")
count_heatmap$no_copper <- log10(rowMeans(count_heatmap[,1:3]))
count_heatmap$copper    <- log10(rowMeans(count_heatmap[,4:6]))

fold_change <- res$log2FoldChange %>%
  as_tibble() %>%
  mutate(genes = rownames(res), log2foldchange = value) %>%
  dplyr::select(genes, log2foldchange)

count_heatmap           <- count_heatmap %>% 
  mutate(genes = rownames(counts(dds_results))) %>%
  dplyr::filter(genes %in% top25) %>%
  dplyr::select(no_copper, copper, genes) 

count_heatmap_mod <- count_heatmap %>%
  gather(genes) %>%
  mutate(geneid = rep(count_heatmap$genes,2)) %>%
  left_join(fold_change, by = c("geneid" = "genes")) %>%
  dplyr::arrange(desc(log2foldchange))

names_list <- list()
id_list <- list()
for(i in 1:length(count_heatmap2$genes)){
  if (length(arcog_pfu$name[arcog_pfu$gene == count_heatmap2$geneid[1]])>0 & length(arcog_pfu$old_name[arcog_pfu$gene == count_heatmap2$geneid[i]])>0){
    names_list[i] <- arcog_pfu$name[arcog_pfu$gene == count_heatmap2$geneid[i]]
    id_list[i] <- arcog_pfu$old_name[arcog_pfu$gene == count_heatmap2$geneid[i]]
  }
  else{
    names_list[i] <- NA
    id_list[i] <- NA
  }
}

#..............................set names
count_heatmap_mod$geneid_extended <- paste(names_list, id_list)

#..............................set colors
colfunc <- colorRampPalette(c("white", "#F1605D"))

#..............................heatmap of log2-transformed counts (Fig. 3d panel 2)
pdf(here("/figures/heatmap_deseq.pdf"), 
    width = 14, height = 7, paper = "special",onefile=FALSE)
ggplot(data = count_heatmap2, aes(y = reorder(geneid_extended, log2foldchange), x = genes, color = value, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#150E37", mid = "white", high = "#F1605D", midpoint = 2) +
  scale_color_gradient2(low = "#150E37", mid = "white", high = "#F1605D", midpoint = 2) +
  theme_Publication_white() +
  ylab("") +
  xlab("") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()

#..............................log2 fold change (Fig. 3d panel 2)
pdf(here("/figures/log2_foldchange.pdf"), 
    width = 14, height = 7, paper = "special",onefile=FALSE)
ggplot(data = count_heatmap2, aes(y = reorder(geneid_extended, log2foldchange), x = "change", color = log2foldchange, fill = log2foldchange)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#F1605D", midpoint = 1.9) +
  scale_color_gradient2(low = "white", high = "#F1605D", midpoint = 1.9) +
  theme_Publication_white() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()

 
#................................ARCOG ANALYSIS
#..............................goseq method/set threshold to 0.05 false discovery rate
fdr.threshold <- 0.05
assayed.genes <- rownames(res)
de.genes <- rownames(res)[ which(res$padj < fdr.threshold & res$log2FoldChange > 1) ]

#..............................get de-genes and background set in a list
gene.vector <-  as.integer(assayed.genes%in%de.genes)
names(gene.vector) <- assayed.genes

#..............................add length information
feature_info <- featureCounts(bam_list[1],verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "exon", allowMultiOverlap = F, isLongRead = F)$annotation %>%
  as_tibble() %>%
  dplyr::filter(GeneID %in% names(gene.vector))

pwf <- nullp(gene.vector, bias.data = feature_info$Length,'ensGene',plot.fit=FALSE)

#..............................add arCOG identifiers
arcog_table_m <- arcog_pfu %>%
  dplyr::filter(gene %in% names(gene.vector)) %>%
  dplyr::select(gene, arCOG) %>%
  as.data.frame()
category.vector <- arcog_table_m$arCOG
names(category.vector) <- as.factor(arcog_table_m$new)

#..............................use goseq to determine differntial arCOG terms
goseq_results <- goseq(pwf, gene2cat = arcog_table_m, use_genes_without_cat=TRUE)

goseq_results_filtered <- goseq_results %>%
  dplyr::filter(over_represented_pvalue < 0.05)

#..............................visualize the results / 1. table of arCOGS with numbers that are differentially regulated after copper shock.

#..............................test set
arcog_table <- arcog_table_m %>%
  group_by(arCOG) %>%
  summarise(counts = n()) %>%
  mutate(set = "test_set")

#..............................background genomic set
arcog_bg <- arcog_pfu %>%
  dplyr::select(gene, arCOG) %>%
  as_tibble() %>%
  group_by(arCOG) %>%
  summarise(counts = n()/1925*25) %>%
  mutate(set = "bg_set")

#..............................join sets
arcogs <- rbind(arcog_table,arcog_bg) 

#..............................radarplot of arCOGs with more than 1 count in category
goseq_results_plotting <- goseq_results %>% 
  dplyr::filter(numDEInCat > 0) %>%
  dplyr::select(numDEInCat, numInCat) %>%
  mutate(numDEInCat_100 = numDEInCat/34,
         numInCat_100 = numInCat/1879) %>%
  dplyr::select(numDEInCat_100, numInCat_100) %>% 
  t() %>%
  as_tibble() 
colnames(goseq_results_plotting) <- goseq_results$category[goseq_results$numDEInCat > 0]
goseq_results_plotting_n <- cbind("names" = c("cu_induced","total"), goseq_results_plotting)

pdf(here("figures/radarchart_arcogs.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggradar(goseq_results_plotting_n, grid.min = 0, grid.max = 0.3, grid.mid = 0.15, label.gridline.max = "20%", gridline.mid.colour = "black", 
        label.gridline.min = "0%", label.gridline.mid = "10%", group.point.size = 4, group.line.width = 0.75, group.colours = c(viridis_pal(option = "magma")(10)[7.5], viridis_pal(option = "magma")(10)[1]))       
dev.off()

