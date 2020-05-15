CopR, a global regulator of transcription to maintain copper homeostasis
in Pyrococcus furiosus
================
Felix Grünberger<sup>1</sup>, Robert Reichelt<sup>1</sup>, Ingrid
Waege<sup>1</sup>, Verena Niedermaier<sup>1</sup>, Korbinian
Bronner<sup>1</sup>, Korbinian Bronner<sup>1</sup>, Marcell
Kaljanac<sup>2</sup>, Nina Weber<sup>1</sup>, Zubeir El
Ahmad<sup>1</sup>, Lena Knauss<sup>1</sup>, Gregor Madej<sup>2</sup>,
Christine Ziegler<sup>2</sup>, Dina Grohmann<sup>1</sup>, Winfried
Hausner<sup>1°</sup>  

<sup>1</sup> Institute of Microbiology and Archaea Centre, University of
Regensburg, Universitätsstraße 31, 93053 Regensburg, Germany

<sup>2</sup> Department of Structural Biology, Institute of Biophysics
and Physical Biochemistry, University of Regensburg, Universitätsstraße
31, 93053 Regensburg, Germany

<sup>°</sup> Corresponding author

  - [Information about this
    repository](#information-about-this-repository)
  - [Data generation](#data-generation)
  - [Data analysis](#data-analysis)
      - [Raw data quality control](#raw-data-quality-control)
      - [Filtering of raw reads](#filtering-of-raw-reads)
  - [Mapping of reads (`STAR`)](#mapping-of-reads-star)
  - [Converting bam files to other
    format](#converting-bam-files-to-other-format)
  - [Preparing count matrices
    (`featurecounts`)](#preparing-count-matrices-featurecounts)
      - [Calculating transcript abundances and differential expression
        analysis](#calculating-transcript-abundances-and-differential-expression-analysis)
          - [Differential gene
            expression](#differential-gene-expression)
          - [ChIP-seq analysis](#chip-seq-analysis)
      - [Data availability](#data-availability)
          - [Raw sequencing files](#raw-sequencing-files)
          - [Additional data](#additional-data)
      - [License](#license)
      - [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

-----

## Information about this repository

This is the repository for the manuscript “CopR, a global regulator of
transcription to maintain copper homeostasis in Pyrococcus furiosus”. It
contains a description of the bioinformatical tools and custom
[Rscripts](Rscripts) to process data and a description of the downstream
analysis.

The repository is currently actively developed.

[![Active
Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)

## Data generation

## Data analysis

### Raw data quality control

Quality control at the raw data level using
<a target="_blank" href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">`FastQC`</a>
was done to examine most common quality parameters (the total number of
reads sequenced, GC content and the overall base quality score). Quality
filtering parameters were based on this analysis.

### Filtering of raw reads

We used
<a href = "http://www.usadellab.org/cms/?page=trimmomatic">`trimmomatic`</a>
(v. 0.36) in single-end-mode to trim raw reads (Bolger, Lohse, and
Usadel ([2014](#ref-Bolger2014))). We start trimming the remaining reads
by certain quality scores (PHRED33 based) with calculations in a sliding
window of 4 bases and a minumum length of 12 bases.

    for file in $data_folder/sortmerna_data/*_rna.fastq
    do 
      xbase=${file##*/}
      filename=${xbase%%.*}
      java -jar trimmomatic-0.36.jar SE -phred33 -threads 4 $file $data_folder"/clean_data/"$filename"_trimmed.fq.gz" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:12
      echo $filename "trimming finished \n\n"
    done

## Mapping of reads (`STAR`)

Reads are mapped to the genome

``` bash
# executable binary in following folder:
cd /Users/felix/Downloads/STAR-2.5.4b/bin/MacOSX_x86_64/
```

``` bash
# generate a index file for the genome based on .gtf file
./STAR.dms \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir $data_folder"/genome_data/star" \
--genomeFastaFiles $data_folder"/genome_data/CP023154.fa" \
--sjdbGTFfile $data_folder"/genome_data/CP023154.gtf" \
--genomeSAindexNbases 10 

# align reads
for file in $data_folder/clean_data/*.fq.gz
do 
  xbase=${file##*/}
  filename=${xbase%%.*}
  mkdir $data_folder"/mapping_data/star/"$filename
  ./STAR \
  --runMode alignReads \
  --runThreadN 4 \
  --readFilesCommand gunzip -c \
  --genomeDir $data_folder"/genome_data/star" \
  --readFilesIn  $file \
  --outFileNamePrefix $data_folder"/mapping_data/star/"$filename/$filename"_" \
  --outSAMtype BAM SortedByCoordinate \
  --limitIObufferSize 300000000 \
  --limitBAMsortRAM 4199333281 
  echo $filename "mapping finished" \n\n 
done
```

## Converting bam files to other format

Convert files to other formats (*BigWig*, *wiggle*) that are usually
used for viewing mapped reads or for further analysis.

Conversion to WIG is done using `bam2wig` script:

``` bash
# in Documents/scripts
for file in $data_folder/mapping_data/star/*/*.bam
do 
  xbase=${file##*/}
  filename=${xbase%%.*}
  ./bam2wig $file
  echo $filename "bam2wig finished" \n\n
done
```

# Preparing count matrices (`featurecounts`)

``` r
# load libraries
library("Rsubread")
library("stringr")
library("data.table")
# load files 
gtf_file <- paste(data_folder,"/genome_data/CP023154.gtf", sep = "")

data_folder <- "/Users/felix/Documents/R/differential_0739/data"
mappedFiles <- dir(file.path(paste(data_folder, "/mapping_data/star/", sep  ="")), pattern="*trimmed", full.name=T)
names(mappedFiles) <- gsub(".bed","",basename(mappedFiles))
mappedFiles


file_folder <- list.files(path = paste(data_folder,"/mapping_data/star", sep = ""), full.names = T)
file_path <- list.files(path = paste(data_folder,"/mapping_data/star", sep = ""), full.names = F)
bam_list <- str_c(file_folder,"/", file_path, "_Aligned.sortedByCoord.out.bam",sep = "") 

# calculate count matrix
counts <- featureCounts(bam_list,verbose = F,annot.ext = gtf_file, isGTFAnnotationFile = T, nthreads = 4,GTF.featureType = "exon", allowMultiOverlap = F, isLongRead = F)$counts
```

After calculating the count matrix, column names are annotated according
to experimental conditions and design.

## Calculating transcript abundances and differential expression analysis

``` r
library("tidyverse")
colnames_data <- c("replicate", "strain", "treatment")
colData <- matrix(ncol = 3, nrow = 12)
colnames(colData) <- colnames_data 
colData <- data.table(colData) %>%
  mutate(replicate = as.factor(rep(1:3,4)),
         strain = as.factor(c(rep(52,6), rep(70,3), rep(74,3))),
         treatment = as.factor(c(rep("no",3), rep("yes",3), rep("no",3), rep("no",3))))
```

Now we construct a *DESeqDataSet* object with:

``` r
library("DESeq2")
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[,c(1:6)], 
                                            colData = colData[c(1:6),],
                                            design = ~treatment)
# prefilter data set: remove rows (genes) with no counts
dds <- dds[ rowSums(counts(dds)) > 1, ]
```

In that case we are comparing strain 52 with 52 copper shock.

### Differential gene expression

### ChIP-seq analysis

## Data availability

### Raw sequencing files

Raw sequence data have been uploaded to the NCBI sequence read archive
(<a href="https://www.ncbi.nlm.nih.gov/sra">SRA</a>) and are available
under project accession number
<a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA603674">PRJNA603674</a>.

### Additional data

-----

## License

This project is under the general MIT License - see the
[LICENSE](LICENSE) file for details

## References

<div id="refs" class="references hanging-indent">

<div id="ref-Bolger2014">

Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A
flexible trimmer for Illumina sequence data.” *Bioinformatics*.
<https://doi.org/10.1093/bioinformatics/btu170>.

</div>

</div>
