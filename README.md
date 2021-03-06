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
  - [Data analysis](#data-analysis)
      - [Raw data quality control](#raw-data-quality-control)
      - [Filtering of raw reads](#filtering-of-raw-reads)
      - [Mapping of reads](#mapping-of-reads)
          - [RNA-seq](#rna-seq)
          - [ChIP-seq](#chip-seq)
      - [Differential gene expression
        analysis](#differential-gene-expression-analysis)
      - [ChIP-seq analysis & integration with RNA-seq
        data](#chip-seq-analysis-integration-with-rna-seq-data)
  - [Data availability](#data-availability)
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

## Data analysis

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
Usadel ([2014](#ref-Bolger2014))).

We start trimming the remaining reads by certain quality scores (PHRED33
based, cut.off Phread score 20) with calculations in a sliding window of
4 bases and a minumum length of 12 bases (40 for ChIP-seq).

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

#### RNA-seq

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

#### ChIP-seq

Reads were mapped to the *P. furiosus* genome using `Bowtie 2` (v.
2.2.3) with default settings (Langmead and Salzberg
[2012](#ref-Langmead2012)).

``` bash
#!/bin/bash

# Indexing a reference genome
cd $data_folder/mapped_data
bowtie2-build $data_folder/genome_data/CP023154.fasta CP023154

# Align reads
for file in $data_folder/trimmed_data/*.fastq.gz
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}  
  (bowtie2 -x CP023154 -U $file -S $data_folder"/mapped_data/"$filename".sam") 2>$data_folder"/mapped_data/"$filename".sam".log
  echo $filename mapping finished
done
```

SAM files were converted to sorted BAM files using samtools and extended
towards the 3´direction to their fragment-size (200 bp size selection
during library preparation) to better represent the precise protein-DNA
interaction (Li et al. [2009](#ref-Li2009)).

``` bash
#!/bin/bash

# sam to indexed sorted bam
for file in $data_folder/mapped_data/*.sam
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  samtools view -bS -q 40 $file > $data_folder"/mapped_data/"$filename."bam"
  samtools sort $data_folder"/mapped_data/"$filename."bam" -o $data_folder"/mapped_data/"$filename."sorted.bam"
  samtools index $data_folder"/mapped_data/"$filename."sorted.bam"
  echo $filename sam_to_bam finished
done

# extend reads towards their 3´end
for file in $data_folder/mapped_data/*.sorted.bam
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  bamToBed -i $file | slopBed -i - -g $data_folder"/genome_data/chrom.sizes.txt" -s -r 150 -l 0 | bedToBam -i - -g $data_folder"/genome_data/chrom.sizes.txt" > $data_folder"/mapped_data/bedtools_extended/"$filename"_extended.bam"
  echo $filename extension finished
done

# sort & index
for file in $data_folder/mapped_data/bedtools_extended/*_extended.bam
do 
  filename_extended=$file##*/
  filename=$filename_extended%%.*
  samtools sort $file -o $data_folder/mapped_data/bam_files/$filename".sorted.extended.bam"
  samtools index $data_folder/mapped_data/bam_files/$filename".sorted.extended.bam"
done

# convert to bedgraph (for log2-calculations in R)
for file in $data_folder/mapped_data/bam_files/*.sorted.extended.bam
do 
  filename_extended=${file##*/}
  filename=${filename_extended%%.*}
  genomeCoverageBed -ibam $file -d > $data_folder"/mapped_data/bed_files/"$filename".sorted.extended.position.bedgraph"
  echo $filename position specific genomecoverage finished
done
```

Position-specific enrichments in BED format were calculated by i)
scaling the reads in each dataset according to sequencing depth, ii)
calculation of the ratio between IP and input for each replicate, and
iii) averaging of the IP/input ratio from the biological triplicates and
taking the log2 for comparison (compare
[`normalise_chipseq`](Rscripts/normalise_chipseq.R)).

### Differential gene expression analysis

We performed differential expression analysis using DESeq2 (Love,
Anders, and Huber [2014](#ref-Love2014)). The analysis performed in R,
including preparation of count matrices, statistical testing and
exploratory data analysis can be retraced using the
[`deseq2_analysis`](Rscripts/deseq2_analysis.R) script.

### ChIP-seq analysis & integration with RNA-seq data

The steps used to perform ChIP-seq enrichment analysis and integration
with RNA-seq data are shown in the
[`chipseq_downstream`](Rscripts/chipseq_downstream.R) script.

## Data availability

Raw sequence data have been uploaded to the NCBI sequence read archive
(<a href="https://www.ncbi.nlm.nih.gov/sra">SRA</a>) and are available
under project accession number
<a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA603674">PRJNA603674</a>.

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

<div id="ref-Langmead2012">

Langmead, Ben, and Steven L. Salzberg. 2012. “Fast gapped-read alignment
with Bowtie 2.” *Nature Methods*. <https://doi.org/10.1038/nmeth.1923>.

</div>

<div id="ref-Li2009">

Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils
Homer, Gabor Marth, Goncalo Abecasis, and Richard Durbin. 2009. “The
Sequence Alignment/Map format and SAMtools.” *Bioinformatics* 25 (16):
2078–9. <https://doi.org/10.1093/bioinformatics/btp352>.

</div>

<div id="ref-Love2014">

Love, M. I., Simon Anders, and Wolfgang Huber. 2014. *Differential
analysis of count data - the DESeq2 package*. Vol. 15. 12.
<https://doi.org/110.1186/s13059-014-0550-8>.

</div>

</div>
