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
<<<<<<< HEAD
install.packages("IHW")
=======

>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ANALYSIS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................DESeq2 pipeline/including featurecounts
<<<<<<< HEAD

#................................genome annotation file
gtf_file <- here("data/genome_data/CP023154.gtf")

#................................mapped files
file_folder <- list.files(path = here("data/mapped_data/star"), full.names = T)
file_path <- list.files(path = here("data/mapped_data/star"), full.names = F)
=======
#................................path
data_folder="/Users/felix/Documents/R/differential_0739/data"

#................................genome annotation file
gtf_file <- paste(data_folder,"/genome_data/CP023154.gtf", sep = "")

#................................mapped files
file_folder <- list.files(path = paste(data_folder,"/mapping_data/star", sep = ""), full.names = T)
file_path <- list.files(path = paste(data_folder,"/mapping_data/star", sep = ""), full.names = F)
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
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
<<<<<<< HEAD
=======
counts_table[counts_table$gene == "PFDSM3638_03680",]
counts_table[counts_table$gene == "PFDSM3638_03675",]
counts_table[counts_table$gene == "PFDSM3638_03690",]
counts_table[counts_table$gene == "PFDSM3638_03585",]

>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0

#................................DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                      colData = colData,
                                      design = ~treatment)

#................................remove genes with no counts
dds <- dds[ rowSums(counts(dds)) > 1, ]

#................................transformation of data / vst (variance stabilizing transformation for negative binomial data with a dispersion-mean trend)
vsd <- vst(dds, blind = T)
<<<<<<< HEAD
results_table <- as_tibble(assay(vsd)) 
colnames(results_table) <- as.factor(1:6) 
table <- results_table %>%  
=======
helper <- as_tibble(assay(vsd)) 
colnames(helper) <- as.factor(1:6) 
table <- helper %>%  
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
  mutate(means_ = rowMeans(.[1:6]),
         sd_ = rowSds(as.matrix(.[1:6])),
         genename = rownames(.))

<<<<<<< HEAD
#................................PCA PLOT (FIG. 3a)
pcaData <- plotPCA(vsd, intgroup = c("replicate", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(here("/figures/pca_plot.pdf"), 
=======
#................................PCA PLOT -- SUPPLEMENTARY FIGURE
pcaData <- plotPCA(vsd, intgroup = c("replicate", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191004_pca_plot_deseq.pdf", 
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
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

<<<<<<< HEAD
#................................BEESWARM PLOT FOR PF0740 (FIG. 3b)
coppertransporter <- rownames(vsd)[which(rownames(vsd) == "PFDSM3638_03695")]
geneCounts <- plotCounts(dds, gene = coppertransporter, intgroup = c("replicate","treatment"), returnData = T)

pdf(here("figures/CopA_counts.pdf"), 
=======
#................................BEESWARM PLOT FOR PF0740 -- MAIN FIGURE
coppertransporter <- rownames(vsd)[which(rownames(vsd) == "PFDSM3638_03695")]
geneCounts <- plotCounts(dds, gene = coppertransporter, intgroup = c("replicate","treatment"), returnData = T)

pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191004_pf0740_counts_deseq.pdf", 
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
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
<<<<<<< HEAD

=======
summary(res)
results(dds_results)
save(res, file = here("data/raw_data/results_deseq2"))
load(here("data/raw_data/results_deseq2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("Hmisc")
library("DESeq2")
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
#................................datapoint with NA of padj have too less counts and are therefore excluded
res_ggplot <- res %>%
  as_tibble() %>%
  rownames_to_column() %>%
  mutate(padj = ifelse(is.na(padj) == T, 1,padj)) %>%
  mutate(type_of_regulation = ifelse(log2FoldChange > 1 & padj < 0.05, "up", 
                                     ifelse(log2FoldChange < -1 & padj < 0.05, "down",
                                            ifelse(is.na(padj), "up","not significant"))),
         id = rownames(res))

<<<<<<< HEAD


#................................MA PLOT OF ALL GENES THAT ARE REGULATED UPON COPPER SHOCK
pdf(here("/figures/ma_plot.pdf"), 
=======
res_up <- res %>%
  as_tibble() %>%
  rownames_to_column() %>%
  mutate(padj = ifelse(is.na(padj) == T, 1,padj),
         id = rownames(res)) %>%
  filter(padj < 0.05, log2FoldChange >= 1) %>%
  left_join(arcog_pfu, by = c("id" = "gene"))

res_down <- res %>%
  as_tibble() %>%
  rownames_to_column() %>%
  mutate(padj = ifelse(is.na(padj) == T, 1,padj),
         id = rownames(res)) %>%
  filter(padj < 0.05, log2FoldChange <= 0) 
arcog_pfu <- fread("/Users/felixgrunberger/Documents/R/arcog_pfu/data/rpm_pyrocoocus_mod_names.tsv")
res_down_n <- res_down %>%
  left_join(arcog_pfu, by = c("id" = "gene"))
res_down_n

sum(res_ggplot$type_of_regulation == "up")
sum(res_ggplot$type_of_regulation == "down")

save(res_ggplot, file = here("/data/raw_data/deseq_results"))
res_ggplot %>%
  dplyr::filter(type_of_regulation == "up", id == "PFDSM3638_03695")

#................................MA PLOT OF ALL GENES THAT ARE REGULATED UPON COPPER SHOCK
pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191004_maplot_deseq.pdf", 
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
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

<<<<<<< HEAD
#................................add arcog information and save to file
arcog_pfu <- fread(here("data/arcog_data/arcog_pfu.tsv"))
res_ggplot_arcog <- res_ggplot %>%
  left_join(arcog_pfu, by = c("id" = "gene"))
write.table(res_ggplot_arcog, file = here("tables/deseq2_table.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")


=======
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
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

<<<<<<< HEAD
count_heatmap_mod <- count_heatmap %>%
=======
count_heatmap2 <- count_heatmap %>%
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
  gather(genes) %>%
  mutate(geneid = rep(count_heatmap$genes,2)) %>%
  left_join(fold_change, by = c("geneid" = "genes")) %>%
  dplyr::arrange(desc(log2foldchange))

<<<<<<< HEAD
=======
#..............................add gene name information
arcog_pfu <- fread("/Users/felixgrunberger/Documents/R/arcog_pfu/data/rpm_pyrocoocus_mod_names.tsv")

>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
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

<<<<<<< HEAD
#..............................set names
count_heatmap_mod$geneid_extended <- paste(names_list, id_list)

#..............................set colors
colfunc <- colorRampPalette(c("white", "#F1605D"))

#..............................heatmap of log2-transformed counts (Fig. 3d panel 2)
pdf(here("/figures/heatmap_deseq.pdf"), 
=======
count_heatmap2$geneid_extended <- paste(names_list, id_list)

colfunc <- colorRampPalette(c("white", "#F1605D"))
colfunc(10)

#..............................log2 counts
pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191004_log2counts_heatmap_deseq.pdf", 
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
    width = 14, height = 7, paper = "special",onefile=FALSE)
ggplot(data = count_heatmap2, aes(y = reorder(geneid_extended, log2foldchange), x = genes, color = value, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#150E37", mid = "white", high = "#F1605D", midpoint = 2) +
  scale_color_gradient2(low = "#150E37", mid = "white", high = "#F1605D", midpoint = 2) +
  theme_Publication_white() +
<<<<<<< HEAD
  ylab("") +
  xlab("") +
=======
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()

<<<<<<< HEAD
#..............................log2 fold change (Fig. 3d panel 2)
pdf(here("/figures/log2_foldchange.pdf"), 
=======
#..............................log2 fold change
pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191004_log2foldchange_heatmap_deseq.pdf", 
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
    width = 14, height = 7, paper = "special",onefile=FALSE)
ggplot(data = count_heatmap2, aes(y = reorder(geneid_extended, log2foldchange), x = "change", color = log2foldchange, fill = log2foldchange)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "#F1605D", midpoint = 1.9) +
  scale_color_gradient2(low = "white", high = "#F1605D", midpoint = 1.9) +
  theme_Publication_white() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0))
dev.off()

<<<<<<< HEAD
 
#................................ARCOG ANALYSIS
=======
log2(105)
log2(290)

#................................ARCOG ANALYSIS
#..............................combine gene names with arcogs in pfu
#arcog <- read_csv("/Volumes/arCOG/ar14.arCOG.csv", skip = 1, col_names = c("pos", "org", "protein_id", "a", "b", "c", "arcog","category", "COG")) %>%
#  dplyr::filter(org == "Pyrococcus_furiosus_DSM_3638_uid57873") %>%
#  dplyr::select(protein_id, arcog,category,COG)
arcog_names <- fread("/Volumes/arCOG/ar14.arCOGdef.tab", skip = 0) %>%
  dplyr::select(V1,V4)
#
#gi2gbk_table <- fread("/Volumes/arCOG/ar14.gi2gbk.tab", skip = 0) 
#
#arcog_combined <- left_join(arcog, arcog_names, by = c("arcog" = "V1")) %>%
#  left_join(gi2gbk_table, by = c("protein_id" = "V1"))
#
## > read in table with pdu to genomic information
#protein_to_name <- fread("/Users/felixgrunberger/Desktop/pfu_dsm3639_arcog_pid.tsv", col.names = c("Strand","Length","PID","Gene","Synonym","Code","COG","Product", "x")) %>%
#  dplyr::select(-x)
#
#arcog_combined_all <- left_join(protein_to_name, arcog_combined, by = c("PID" = "protein_id")) %>%
#  dplyr::select(PID, Synonym, Product, V4,arcog)
#
## >  new names?
#old_new <- read_excel("/Users/felixgrunberger/Downloads/Table_2_Next Generation DNA-Seq and Differential RNA-Seq Allow Re-annotation of the Pyrococcus furiosus DSM 3638 Genome and Provide Insights Into Arch.XLSX", skip = 3, col_names = c("new", "start", "end", "old_r", "old_rs", "old_rse", "old", "new2", "name", "x")) %>%
#  dplyr::select(new, start, end, old, name)
#
#arcog_all <- left_join(arcog_combined, protein_to_name, by = c("protein_id" = "PID")) %>%
#  left_join(old_new, by = c("Synonym" = "old")) 
#save(arcog_all, file = here("data/raw_data/arcog_pfu_table"))
load(here("data/raw_data/arcog_pfu_table"))

>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
#..............................goseq method/set threshold to 0.05 false discovery rate
fdr.threshold <- 0.05
assayed.genes <- rownames(res)
de.genes <- rownames(res)[ which(res$padj < fdr.threshold & res$log2FoldChange > 1) ]

#..............................get de-genes and background set in a list
<<<<<<< HEAD
gene.vector <-  as.integer(assayed.genes%in%de.genes)
names(gene.vector) <- assayed.genes

#..............................add length information
feature_info <- featureCounts(bam_list[1],verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "exon", allowMultiOverlap = F, isLongRead = F)$annotation %>%
  as_tibble() %>%
  dplyr::filter(GeneID %in% names(gene.vector))
=======
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

#..............................add length information
annotation
feature_info <- featureCounts(bam_list[1],verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "exon", allowMultiOverlap = F, isLongRead = F)$annotation %>%
  as_tibble() 
dplyr::filter(GeneID %in% names(gene.vector))
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0

pwf <- nullp(gene.vector, bias.data = feature_info$Length,'ensGene',plot.fit=FALSE)

#..............................add arCOG identifiers
<<<<<<< HEAD
arcog_table_m <- arcog_pfu %>%
  dplyr::filter(gene %in% names(gene.vector)) %>%
  dplyr::select(gene, arCOG) %>%
  as.data.frame()
category.vector <- arcog_table_m$arCOG
names(category.vector) <- as.factor(arcog_table_m$new)

#..............................use goseq to determine differntial arCOG terms
goseq_results <- goseq(pwf, gene2cat = arcog_table_m, use_genes_without_cat=TRUE)
=======
category_mapping <- arcog_all %>%
  dplyr::rename(arCOG = COG.x) %>%
  dplyr::filter(new %in% names(gene.vector)) %>%
  dplyr::select(new, arCOG) %>%
  as.data.frame()
category.vector <- category_mapping$arCOG
names(category.vector) <- as.factor(category_mapping$new)

#..............................use goseq to determine differntial arCOG terms
goseq_results <- goseq(pwf, gene2cat = category_mapping, use_genes_without_cat=TRUE)
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0

goseq_results_filtered <- goseq_results %>%
  dplyr::filter(over_represented_pvalue < 0.05)

#..............................visualize the results / 1. table of arCOGS with numbers that are differentially regulated after copper shock.
<<<<<<< HEAD

#..............................test set
arcog_table <- arcog_table_m %>%
=======
#. 2 categories and spiderplot

#..............................filter out 25 most significantly regulated genes
arcog_table <- arcog_all %>%
  dplyr::rename(arCOG = COG.x) %>%
  dplyr::filter(new %in% rownames(assay(vsd)[ top25, ])) %>%
  dplyr::select(arCOG) %>%
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
  group_by(arCOG) %>%
  summarise(counts = n()) %>%
  mutate(set = "test_set")

#..............................background genomic set
<<<<<<< HEAD
arcog_bg <- arcog_pfu %>%
  dplyr::select(gene, arCOG) %>%
  as_tibble() %>%
=======
arcog_bg <- arcog_all %>%
  dplyr::select(new, COG.x) %>%
  dplyr::rename(arCOG = COG.x) %>%
  as.tibble() %>%
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
  group_by(arCOG) %>%
  summarise(counts = n()/1925*25) %>%
  mutate(set = "bg_set")

#..............................join sets
arcogs <- rbind(arcog_table,arcog_bg) 

<<<<<<< HEAD
=======
#..............................plot
ggplot(arcogs, aes(x = arCOG, y = as.numeric(counts), fill = set)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.5), alpha = 0.95) +
  coord_flip() +
  theme_bw() +
  scale_fill_viridis(discrete =T,option = "magma", begin = 0.7, end = 0.1)

sum(goseq_results$numDEInCat)
sum(goseq_results$numInCat)
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
#..............................radarplot of arCOGs with more than 1 count in category
goseq_results_plotting <- goseq_results %>% 
  dplyr::filter(numDEInCat > 0) %>%
  dplyr::select(numDEInCat, numInCat) %>%
  mutate(numDEInCat_100 = numDEInCat/34,
         numInCat_100 = numInCat/1879) %>%
  dplyr::select(numDEInCat_100, numInCat_100) %>% 
  t() %>%
<<<<<<< HEAD
  as_tibble() 
colnames(goseq_results_plotting) <- goseq_results$category[goseq_results$numDEInCat > 0]
goseq_results_plotting_n <- cbind("names" = c("cu_induced","total"), goseq_results_plotting)

pdf(here("figures/radarchart_arcogs.pdf"), 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggradar(goseq_results_plotting_n, grid.min = 0, grid.max = 0.3, grid.mid = 0.15, label.gridline.max = "20%", gridline.mid.colour = "black", 
        label.gridline.min = "0%", label.gridline.mid = "10%", group.point.size = 4, group.line.width = 0.75, group.colours = c(viridis_pal(option = "magma")(10)[7.5], viridis_pal(option = "magma")(10)[1]))       
=======
  as_tibble() %>%
  mutate(O = c(0.14102564, 0.02927089))
colnames(goseq_results_plotting)[1:10] <- goseq_results$category[goseq_results$numDEInCat > 0]

pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191002_radarchart_arcogs.pdf", 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggradar(goseq_results_plotting, grid.min = 0, grid.max = 0.3, grid.mid = 0.15, label.gridline.max = "20%", gridline.mid.colour = "black", 
        label.gridline.min = "0%", label.gridline.mid = "10%", group.point.size = 4, group.line.width = 0.75, group.colours = c(viridis_pal(option = "magma")(10)[1], viridis_pal(option = "magma")(10)[7.5]))       
dev.off()

#### SAME FOR DOWNREGULATED GENES
res_down_n
#..............................goseq method/set threshold to 0.05 false discovery rate
fdr.threshold <- 0.05
assayed.genes <- rownames(res)
de.genes <- rownames(res)[ which(res$padj < fdr.threshold & res$log2FoldChange <= 0) ]

#..............................get de-genes and background set in a list
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes

#..............................add length information
feature_info <- featureCounts(bam_list[1:6],verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "exon", allowMultiOverlap = F, isLongRead = F)$annotation %>%
  as_tibble() %>%
  dplyr::filter(GeneID %in% names(gene.vector))

pwf <- nullp(gene.vector, bias.data = feature_info$Length,'ensGene',plot.fit=FALSE)

#..............................add arCOG identifiers
category_mapping <- arcog_all %>%
  dplyr::rename(arCOG = COG.x) %>%
  dplyr::filter(new %in% names(gene.vector)) %>%
  dplyr::select(new, arCOG) %>%
  as.data.frame()
category.vector <- category_mapping$arCOG
names(category.vector) <- as.factor(category_mapping$new)

#..............................use goseq to determine differntial arCOG terms
goseq_results <- goseq(pwf, gene2cat = category_mapping, use_genes_without_cat=TRUE)

goseq_results_filtered <- goseq_results %>%
  dplyr::filter(over_represented_pvalue < 0.05)

#..............................radarplot of arCOGs with more than 1 count in category

goseq_results_plotting <- goseq_results %>% 
  dplyr::filter(numDEInCat > 0) %>%
  dplyr::select(numDEInCat, numInCat) %>%
  mutate(numDEInCat_100 = numDEInCat/sum(goseq_results$numDEInCat),
         numInCat_100 = numInCat/sum(goseq_results$numInCat)) %>%
  dplyr::select(numDEInCat_100, numInCat_100) %>% 
  t() %>%
  as_tibble() %>%
  mutate(O = c(0.14102564, 0.02927089))
colnames(goseq_results_plotting)[1:10] <- goseq_results$category[goseq_results$numDEInCat > 0]

pdf("/Users/felixgrunberger/Documents/R/publication_CopR/figures/raw_plots/191002_radarchart_arcogs.pdf", 
    width = 7, height = 7, paper = "special",onefile=FALSE)
ggradar(goseq_results_plotting, grid.min = 0, grid.max = 0.3, grid.mid = 0.15, label.gridline.max = "20%", gridline.mid.colour = "black", 
        label.gridline.min = "0%", label.gridline.mid = "10%", group.point.size = 4, group.line.width = 0.75, group.colours = c(viridis_pal(option = "magma")(10)[1], viridis_pal(option = "magma")(10)[7.5]))       
>>>>>>> bf71ccd034a73be8be45a21620a7bc51a252f5b0
dev.off()

