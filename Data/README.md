# Gene Specific FLAME's Toy data

These datasets are meant for testing the pipeline, testing parameters and for educational purposes. 

The toy data directory only contains:
- The filtered BED12 reads from the C666-1 long-read sequencing (Decompress before use).
- The standard EBV (NC_007605.1) GTF reference.

These datasets are the bare necessity for testing the ***Gene-Specific*** FLAME module. Toy data for FLAME-GLOW is planned to be amended in a future patch.

Additional compliments for further in-depth analysis using short-read confirmation and canonical splice junction detection needs to be downloaded from public repositories. 

These repositories are:

- NCBI's [EBV (NC_007605.1) reference assembly sequence](https://www.ncbi.nlm.nih.gov/nuccore/NC_007605.1?report=fasta) in fasta-format which are publically available.

- [Edwards, R.H., et al]( https://doi.org/10.1371/journal.ppat.1008071)'s, short-read sequencing data. 

  This data is publically available under the NCBI BioProject ID [PRJNA501807](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA501807). 
  
  The dataset with SRA Accession ID:[SRR8133536](https://www.ncbi.nlm.nih.gov/sra/SRX4954558) was used and is available to download through NCBI's [SRA-tools](https://github.com/ncbi/sra-tools) toolkit.

## Example on running Gene Specific FLAME with the Toy data
This example will use the EBV non-coding RNA RPMS1 as the target gene.
- Using Gene Specific FLAME with minimal addition will use the following syntax:  
 ```sh
 ./FLAME.py -I C666_LR_01.bed -GTF Epstein-BarrVirus.gtf -G RPMS1
 ```
 
- Using Gene Specific FLAME with canonical splice junction detection will use the following syntax:  
 ```sh
 ./FLAME.py -I C666_LR_01.bed -GTF Epstein-BarrVirus.gtf -G RPMS1 -R [Reference.fasta]
 ```
 
- Using Gene Specific FLAME with shortread confirmation will use the following syntax:  
 ```sh
 ./FLAME.py -I C666_LR_01.bed -GTF Epstein-BarrVirus.gtf -G RPMS1 -B [Shortread.bam]
 ```
## Example on running FLAME-GLOW with the Toy data
Currently no toy data exists for the FLAME-GLOW pipeline. <br>
This is planned to be amended in a future patch.


# FLAME-GLOW's Toy data
Currently, there are no toy data to test the FLAME-GLOW module uploaded. <br>
This is planned to be amended in a future patch.
