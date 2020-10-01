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
  
  left_join(vroom(binnedFiles[i],col_names = F), vroom(binnedFiles[i+1],col_names = F), by = c("X1", "X2")) %>%
    left_join(vroom(binnedFiles[i+2],col_names = F), by = c("X1", "X2")) %>%
    left_join(vroom(binnedFiles[i+3],col_names = F), by = c("X1", "X2")) %>%
    left_join(vroom(binnedFiles[i+4],col_names = F), by = c("X1", "X2")) %>%
    left_join(vroom(binnedFiles[i+5],col_names = F), by = c("X1", "X2")) %>%
    mutate(X3.x = .[[3]]/sum(.[[3]])*1000000,
           X3.y = .[[4]]/sum(.[[4]])*1000000,
           X3.x.x = .[[5]]/sum(.[[5]])*1000000,
           X3.y.y = .[[6]]/sum(.[[6]])*1000000,
           X3.x.x.x = .[[7]]/sum(.[[7]])*1000000,
           X3.y.y.y = .[[8]]/sum(.[[8]])*1000000) %>%
    mutate(ip_input_1 = X3.y.y/X3.x,
           ip_input_2 = X3.x.x.x/X3.y,
           ip_input_3 = X3.y.y.y/X3.x.x) %>%
    mutate(ip_input_mean = rowMeans(.[,9:11]),
           counts_log2 = log2(ip_input_mean)) %>%
    dplyr::select(X1, X2, counts_log2)
}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# READ IN DATA & NORMALISE
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# path to extended bed files
exbedFiles_list <- here("/data/mapped_data/bed_files")
exbedFiles <- list.files(path = paste(exbedFiles_list), pattern="*.sorted.extended.position.bedgraph$", full.name=T)
names(exbedFiles) <- gsub(".bedgraph","",basename(exbedFiles))

# calc log2 IP vs input normal
chip739_exbed <- ChIP_norm(chip_project = "739") %>%
  as.tibble() %>%
  mutate(col2 = X2-1,
         col3 = X2,
         col4 = counts_log2) %>%
  dplyr::select(X1, col2, col3, col4)

# calc log2 IP vs input copper-treated
copper_chip739_exbed <- ChIP_norm(chip_project = "copper_739") %>%
  as.tibble() %>%
  mutate(col2 = X2-1,
         col3 = X2,
         col4 = counts_log2) %>%
  dplyr::select(X1, col2, col3, col4)

# export log2-files
write.table(chip739_exbed,quote = F, sep = "\t", row.names = F, col.names = F,
            file = paste(exbedFiles_list, "/log2_739.new.bedgraph", sep = ""))

write.table(copper_chip739_exbed,quote = F, sep = "\t", row.names = F, col.names = F,
            file = paste(exbedFiles_list, "/log2_739_copper.new.bedgraph", sep = ""))

