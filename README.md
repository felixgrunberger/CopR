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
      - [Navigation](#navigation)
      - [Raw data quality control](#raw-data-quality-control)
      - [Filtering of raw reads](#filtering-of-raw-reads)
      - [Mapping of reads](#mapping-of-reads)
      - [Diffential gene expression
        analysis](#diffential-gene-expression-analysis)
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

### Navigation

To follow the data analysis we managed the input/output folders in the
following way:

``` bash
CopR/
├── data/
|   ├── genome_data
|   ├── raw_data
|   ├── clean_data
|   ├── mapped_data
|   ├── fastq_data
|   ├── meme_data
|   └── operon_data
├── Rscrips
├── figures
├── tables/
|   ├── tss_tables
|   ├── tts_tables
|   ├── tu_tables
|   └── counts_tables
├── LICENSE
└── README
```

### Raw data quality control

Quality control at the raw data level using
<a target="_blank" href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">`FastQC`</a>
was done to examine most common quality parameters (total number of
reads sequenced, GC content and the overall base quality score). Quality
filtering parameters were based on this analysis.

### Filtering of raw reads

We used
<a href = "http://www.usadellab.org/cms/?page=trimmomatic">`trimmomatic`</a>
(v. 0.36) in single-end-mode to trim raw reads (Bolger, Lohse, and
Usadel ([2014](#ref-Bolger2014))). We start trimming the remaining reads
by certain quality scores (PHRED33 based) with calculations in a sliding
window of 4 bases and a minumum length of 12 bases.

``` bash
#!/bin/bash

for file in $data_folder/raw_data/*_rna.fastq
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  java -jar trimmomatic-0.36.jar SE -phred33 -threads 4 $file $data_folder"/clean_data/"$filename"_trimmed.fq.gz" LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:12
  echo $filename "trimming finished"
done
```

### Mapping of reads

After trimming, reads were mapped to the genome using `STAR` (v.2.5.4)
(Dobin et al. [2013](#ref-Dobin2013)).

``` bash
#!/bin/bash

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
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  mkdir $data_folder"/mapped_data/star/"$filename
  ./STAR \
  --runMode alignReads \
  --runThreadN 4 \
  --readFilesCommand gunzip -c \
  --genomeDir $data_folder"/genome_data/star" \
  --readFilesIn  $file \
  --outFileNamePrefix $data_folder"/mapped_data/star/"$filename/$filename"_" \
  --outSAMtype BAM SortedByCoordinate \
  --limitIObufferSize 300000000 \
  --limitBAMsortRAM 4199333281 
  echo $filename "mapping finished" 
done
```

After mapping files were converted to other formats (*BigWig*, *wiggle*)
that can be used used for viewing in a genome browser or for further
analysis.  
Strand-specific wig and bigwig files were finally created using
`bam2wig` (Version 1.5, <https://github.com/MikeAxtell/bam2wig>).

``` bash
#!/bin/bash

for file in $data_folder/mapped_data/star/*/*.bam
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  ./bam2wig $file
  echo $filename "bam2wig finished" 
done
```

### Diffential gene expression analysis

We performed differential expression analysis using DESeq2 (Love,
Anders, and Huber [2014](#ref-Love2014)). The analysis performed in R,
including preparation of count matrices, statistical testing and
exploratory data analysis can be retracted using the
[`deseq2_analysis`](Rscripts/deseq2_analysis.R) script.

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

<div id="ref-Dobin2013">

Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow,
Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R.
Gingeras. 2013. “STAR: Ultrafast universal RNA-seq aligner.”
*Bioinformatics* 29 (1): 15–21.
<https://doi.org/10.1093/bioinformatics/bts635>.

</div>

<div id="ref-Love2014">

Love, M. I., Simon Anders, and Wolfgang Huber. 2014. *Differential
analysis of count data - the DESeq2 package*. Vol. 15. 12.
<https://doi.org/110.1186/s13059-014-0550-8>.

</div>

</div>
