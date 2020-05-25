###########################################################################
###########################################################################
###
### LOAD EXTENDED ChIP-SEQ READ FILES AND CALCULATE LOG2-ENRICHMENT
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

#...................................normalize chipseq data
ChIP_norm <- function(chip_project = c("739", "copper_739")){
  binnedFiles <- exbedFiles
  if(chip_project == "739"){
    i <- 1
    
    
  } else if(chip_project == "copper_739"){
    i <- 1+6
  } 
  
  left_join(fread(binnedFiles[i]), fread(binnedFiles[i+1]), by = c("V1", "V2")) %>%
    left_join(fread(binnedFiles[i+2]), by = c("V1", "V2")) %>%
    left_join(fread(binnedFiles[i+3]), by = c("V1", "V2")) %>%
    left_join(fread(binnedFiles[i+4]), by = c("V1", "V2")) %>%
    left_join(fread(binnedFiles[i+5]), by = c("V1", "V2")) %>%
    mutate(V3.x = .[[3]]/sum(.[[3]])*1000000,
           V3.y = .[[4]]/sum(.[[4]])*1000000,
           V3.x.x = .[[5]]/sum(.[[5]])*1000000,
           V3.y.y = .[[6]]/sum(.[[6]])*1000000,
           V3.x.x.x = .[[7]]/sum(.[[7]])*1000000,
           V3.y.y.y = .[[8]]/sum(.[[8]])*1000000) %>%
    mutate(counts_input = rowMeans(.[,3:5]),
           counts_chip = rowMeans(.[,c(6,8)]),
           counts_log2 = log2(counts_chip/counts_input)) %>%
    dplyr::select(V1, V2, counts_log2)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# READ IN DATA & NORMALISE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# path to extended bed files
exbedFiles_list <- here("/data/mapped_data/bed_files")
exbedFiles <- list.files(path = paste(exbedFiles_list), pattern="*.sorted.new.extended.position.bedgraph$", full.name=T)
names(exbedFiles) <- gsub(".bedgraph","",basename(exbedFiles))

# calc log2 IP vs input normal
chip739_exbed <- ChIP_norm(chip_project = "739") %>%
  as.tibble() %>%
  mutate(col2 = V2-1,
         col3 = V2,
         col4 = counts_log2) %>%
  dplyr::select(V1, col2, col3, col4)

# calc log2 IP vs input copper-treated
copper_chip739_exbed <- ChIP_norm(chip_project = "copper_739") %>%
  as.tibble() %>%
  mutate(col2 = V2-1,
         col3 = V2,
         col4 = counts_log2) %>%
  dplyr::select(V1, col2, col3, col4)

# export log2-files
write.table(chip739_exbed,quote = F, sep = "\t", row.names = F, col.names = F,
            file = paste(exbedFiles_list, "/log2_739.bedgraph", sep = ""))

write.table(copper_chip739_exbed,quote = F, sep = "\t", row.names = F, col.names = F,
            file = paste(exbedFiles_list, "/log2_739_copper.bedgraph", sep = ""))

