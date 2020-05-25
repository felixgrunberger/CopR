###########################################################################
###########################################################################
###
### DOWNSTREAM CHIPSEQ ANALYSIS
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

#...................................subsampling intergenic regions
sampleString = function(string) {
  nStart = sample(1:(nchar(string) - 61),1)
  substr(string, nStart, nStart + 60)
}

#...................................make consensus matrix
make_matrix_c <- function(input_term_sequence){
  consensusMatrix(input_term_sequence, as.prob = T) %>%
    t() %>%
    as_tibble() %>%
    mutate(position = c(-60:0)) %>%
    gather(key = position) %>%
    dplyr::rename(base = 1) %>%
    mutate(position = rep(c(-60:0),4))
}

#...................................make read count matrix
make_read_matrix <- function(left_border, right_border, dataset_location){
  dataset <- fread(dataset_location) %>%
    dplyr::select(V4)
  line <- seq(left_border, right_border, by = 1)
  # plus strand
  matrix_p <- matrix(nrow = length(p_plus$V1), ncol = length(line),data = NA)
  for(i in 1:length(p_plus$V1)){
    matrix_p[i,] <- unlist(dataset[p_plus$V2[i] + line])
  }
  # minus strand
  matrix_m <- matrix(nrow = length(p_minus$V1), ncol = length(line),data = NA)
  for(i in 1:length(p_minus$V1)){
    matrix_m[i,] <- unlist(dataset[p_minus$V2[i] - line])
  }
  
  matrix <- rbind(matrix_p, matrix_m[,c(ncol(matrix_m):1)])
  matrix[!is.finite(matrix)] <- NA
  matrix
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ANALYSIS - ENRICHMENT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................location to bed files
data_folder <- here("/data/mapped_data/bed_files")
mappedFiles <- dir(file.path(paste(data_folder, "/mapped_data/bed_files", sep  ="")), pattern="*.sorted.extended.position.bedgraph$", full.name=T)

#...................................merge normal with copper conditions
allmappedFiles <- left_join(fread(mappedFiles[2]), fread(mappedFiles[1]), by = c("V1", "V2", "V3")) 

#...................................genomic coordinates of top10 genes found by differential RNA seq
res <- fread(here("tables/deseq2_table.tsv"))
selected_genes <- res$id[order(res$log2FoldChange, decreasing = T)[1:14]]

gff_table <- readGFF(here("data/genome_data/CP023154.gff")) %>%
  as.tibble() %>%
  dplyr::filter(type == "mRNA") %>%
  mutate(ID = str_sub(ID, 1, 15)) %>%
  dplyr::select(ID, start, end, strand) %>%
  dplyr::filter(ID %in% selected_genes)

x <- min(gff_table$start) - 10000
y <- max(gff_table$end) + 10000

#...................................prepare plot
line <- c(x:y)
chip_line         <- allmappedFiles$V4.x[x:y]
copper_chip_line  <- allmappedFiles$V4.y[x:y]
data_TLE <- matrix(ncol = 3, nrow = y+1-x)
data_TLE[,1] <- line
data_TLE[,2] <- chip_line  
data_TLE[,3] <- copper_chip_line  
data_TLE <- as.data.table(data_TLE)
colnames(data_TLE) <- c("position", "chip", "copper_chip")

#...............................plot with annotation of top10 upregulated genes (Fig. 4a)
pdf(here("/figures/chipseq_genomicregion_top10.pdf"), 
    width = 14, height = 5, paper = "special",onefile=FALSE)
ggplot(data = data_TLE) + 
  geom_line(aes(x = as.numeric(position), y= as.numeric(chip)), size = 1, alpha = 1,
            color = viridis_pal(option="magma")(10)[1]) +
  geom_line(aes(x = as.numeric(position), y= as.numeric(copper_chip)), size = 1, alpha = 1,
            color = viridis_pal(option="magma")(10)[7]) +
  geom_hline(size = 1.5, color = "black", yintercept = 0, linetype = "dashed") +
  theme_bw() + 
  labs(x = "", y = "ChIP occupancy, log2(IP/input)")  +
  scale_x_continuous(expand = c(0,0)) +
  theme(text=element_text(size=30),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rect(data = gff_cds, aes(xmin = start, xmax = end, ymin = -3.5, ymax = -2.5, color = group, fill = group),alpha = 0.5) +
  scale_color_manual(values = c("grey60", "black")) +
  scale_fill_manual(values = c("grey60", "black")) 
dev.off()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA INTEGRATION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................read in fasta
pfu_fasta <- readDNAStringSet(here("/data/genome_data/CP023154.fasta"))
names(pfu_fasta) <- "CP023154"

#...................................transform to data.table and split strands
topgenes <- res %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble

gff_table <- readGFF(here("/data/genome_data/CP023154.gff")) %>%
  as.tibble() %>%
  dplyr::filter(type == "mRNA") %>%
  mutate(ID = str_sub(ID, 1, 15)) %>%
  dplyr::select(ID, start, end, strand)

gff_p <- left_join(topgenes, gff_table, by = c("rowname" = "ID")) %>%
  dplyr::filter(strand == "+") %>%
  mutate(chr = "CP023154", 
         name = "CDS") %>% 
  dplyr::select(chr, start, end, name)

gff_m <- left_join(topgenes, gff_table, by = c("rowname" = "ID")) %>%
  dplyr::filter(strand == "-") %>%
  mutate(chr = "CP023154", 
         start = end,
         name = "CDS") %>% 
  dplyr::select(chr, start, end, name)

#...................................write CDS start position to table 
write.table(gff_p,  file = paste(here(),"/data/raw_data/CP023154_CDS_filtered_plus.gff.bed",sep = ""),row.names = F, col.names = F, quote = F, sep = "\t")     
write.table(gff_m,  file = paste(here(),"/data/raw_data/CP023154_CDS_filtered_minus.gff.bed",sep = ""),row.names = F, col.names = F, quote = F, sep = "\t")     


#...................................read in CDS files
p_plus  <- fread(input = paste(here(),"/data/raw_data/CP023154_CDS_filtered_plus.gff.bed", sep = ""))     
p_minus <- fread(input = paste(here(),"/data/raw_data/CP023154_CDS_filtered_minus.gff.bed", sep = ""))

# --> remove first value
p_plus <- p_plus[-1,]

#...................................set left and right border
x <- 300

#...................................make read matrix from extended bedgraph files
matrix739           <- make_read_matrix(left_border = -x, right_border = x, dataset_location = "/Volumes/TOSHIB/backup_lacie_190320/chip_seq/20181113/mapped_data/bed_files/log2_739_individual_extended.bedgraph")
matrix739_copper    <- make_read_matrix(left_border = -x, right_border = x, dataset_location = "/Volumes/TOSHIB/backup_lacie_190320/chip_seq/20181113/mapped_data/bed_files/log2_739_copper_individual_extended.bedgraph")

#...................................filter genes for strand
topgenes_gff_p_list <- left_join(topgenes, gff_table, by = c("rowname" = "ID")) %>%
  dplyr::filter(strand == "+") %>%
  mutate(chr = "CP023154") 
topgenes_gff_p_list <- topgenes_gff_p_list[-1,]
topgenes_gff_m_list <- left_join(topgenes, gff_table, by = c("rowname" = "ID")) %>%
  dplyr::filter(strand == "-") %>%
  mutate(chr = "CP023154")
topgenes_list <- rbind(topgenes_gff_p_list, topgenes_gff_m_list)

#...................................add count values to loaded rnaseq table
occupancy_table <- matrix(ncol = 3, nrow = length(res$rowname))
occupancy_table[,1] <- topgenes_list$rowname
occupancy_table[,2] <- as.numeric(rowMeans(matrix739))
occupancy_table[,3] <- as.numeric(rowMeans(matrix739_copper))
occupancy_table <- as_tibble(occupancy_table)
colnames(occupancy_table) <- c("genename", "chip739_occupancy", "chip739_copper_occupancy")

#...................................explicit grouping/ gene has to be enriched in both chip seq sets, a significant upregulation & 739 bound under both conditions
chip_and_rnaseq_table <- occupancy_table %>%
  mutate(chip739_occupancy = as.numeric(chip739_occupancy),
         chip739_copper_occupancy = as.numeric(chip739_copper_occupancy)) %>%
  inner_join(res_ggplot, by = c("genename" = "id")) %>%
  mutate(group = ifelse(log2FoldChange > 1 & padj < 0.05 & chip739_occupancy > 1 & chip739_copper_occupancy > 1, "up", "else")) %>%
  left_join(arcog_pfu, by = c("genename" = "gene"))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MOTIF ANALYSIS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
x <- 300
chip_and_rnaseq_table_seq <- left_join(chip_and_rnaseq_table %>% dplyr::filter(group == "up"), gff_table, by = c("genename" = "ID")) %>%
  rowwise() %>%
  mutate(sequence = ifelse(strand == "+", as.character(pfu_fasta$CP023154[(start-x):(start+0)]),
                           as.character(reverseComplement(pfu_fasta$CP023154[(end-0):(end+x)]))))

#...................................write sequences to fasta file for meme analysis
write.fasta(as.list(chip_and_rnaseq_table_seq$sequence), 
            as.list(chip_and_rnaseq_table_seq$old_name), file = here("data/meme_data/chip739_enriched_rna_and_chip.fasta")) 

#...................................MEME ANALYSIS (use intergenic background file derived from all intergenic sequences)
# meme /data/meme_data/chip739_enriched_rna_and_chip.fasta -dna -oc data/meme_data/chip739_enriched_rna_and_chip_intergenicbg -nostatus -time 18000 -mod zoops -nmotifs 5 -minw 3 -maxw 30 -maxsize 200000 -bfile data/genome_data/pfu_intergenic_bg


#...................................MOTIF PLOTTING
color_scale = make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                              cols=viridis_pal(option = "magma")(10)[c(7,1,2,9)])
motif <- read.table(here("data/meme_data/copr_motif_meme.txt"), fill = T, quote = "", sep = "\t") %>%
  dplyr::filter(!grepl("offset", V1)) %>%
  mutate(V1 = as.character(V1))

pdf(here("figures/meme_motif.pdf"), 
    width = 7, height = 3, paper = "special",onefile=FALSE)
ggplot() + 
  geom_logo(motif, font = "helvetica_bold", col_scheme = color_scale, seq_type = "dna") +
  theme_logo() +
  theme_Publication_white() +
  theme(panel.grid.major = element_line(colour = NA),
        axis.ticks.x = element_line(colour = NA), 
        axis.text.x = element_text(size = 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,2), expand = c(0,0))
dev.off()



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# NUCLEOTIDE ENRICHMENT ANALYSIS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................set manually corrected TSS
chip_and_rnaseq_table_seq$custom_TSS <- c(702687,
                                          703398,
                                          704525,
                                          705149,
                                          NA,
                                          705495,
                                          706012,
                                          706274,
                                          716385,
                                          718299,
                                          692601,
                                          NA,
                                          693206,
                                          702602,
                                          716262)
x <- 60
x <- 100
chip_and_rnaseq_table_seq_filtered <- chip_and_rnaseq_table_seq %>%
  dplyr::filter(!is.na(custom_TSS)) %>%
  rowwise() %>%
  mutate(sequence = ifelse(strand == "+", as.character(pfu_fasta$CP023154[(custom_TSS-x):(custom_TSS+0)]),
                           as.character(reverseComplement(pfu_fasta$CP023154[(custom_TSS-0):(custom_TSS+x)]))))

pfu_terminator_matrix <- make_matrix_c(chip_and_rnaseq_table_seq_filtered$sequence)
pfu_background_matrix <- make_matrix_c(pfu_tss$sequence) #tss from pyrococcus with good promoter
pfu_compare_matrix    <- pfu_terminator_matrix %>%
  dplyr::rename(value_terminator = 2) %>%
  left_join(pfu_background_matrix) %>%
  mutate(value_terminator = value_terminator,
         value = value) %>%
  mutate(value = (value_terminator/value))

pdf(here("figures/nucleotide_enrichment_tss_motif_pfu_filtered.pdf"),
    width = 14, height = 7, paper = "special",onefile=FALSE)
ggplot(data = pfu_compare_matrix, aes(x = position, y = (value), color = base)) +
  geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed") +
  geom_smooth(span = 0.09, se = F, size = 2) +
  xlab("Position relative to TSS") +
  ylab("Nucleotide enrichment \n(log2-fold enrichment)") +
  theme_Publication_white() +
  scale_color_viridis_d(option = "magma", begin = 0.8, end = 0) +
  scale_x_continuous(expand = c(0,0))
dev.off()

