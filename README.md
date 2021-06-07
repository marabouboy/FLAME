# FLAME: Full-Length Adjacency Matrix Enumeration.
FLAME:fire: is a long-read splice variant annotation tool.

FLAME has two modules:
- Gene Specific FLAME:  
  - Filtering, translation, quantification of splice variants and exon connectivity through an [Adjacency Matrix](https://en.wikipedia.org/wiki/Adjacency_matrix)
  - Detection of novel exons sites or exon variants through frequency
  - Confirmation of said novel exons sites through the use of:
    - Detection of possible adjacent splice site signals
    - Short-read sequencing reads
- FLAME-GLOW:
  -  Global gene iteration
  -  Filtering and flagging of genes with incongruent annotations
  -  Translation and quantification of genes with sufficient reference annotation

## Requirements
- Python 3.6
- pysam

## Installation

OS X & Linux:

```
git clone https://github.com/marabouboy/FLAME
```

## FLAME usage example
- Gene Specific FLAME under minimal condition adheres to the following syntax:  
  ```sh
  ./FLAME.py -I [INPUT.bed] -GTF [Annotation.gtf] -G [Gene]
  ```
  Minimal Gene Specific FLAME will run:  
  1. *create.ref*  
  2. *filter*  
  3. *translate*  
  4. *quantify*  
  5. *annotated.adjmtx*
  6. *frequency.site*
  7. *frequency.thresh*
  8. *incongruent.adjmtx*  
  
- Gene Specific FLAME with maximal condition adheres to the following syntax:  
  ```sh
  ./FLAME.py -I [INPUT.bed] -GTF [Annotation.gtf] -G [Gene] -R [Reference.fasta] -B [Shortread.bam]
  ```  
  Maximal Gene Specific FLAME will run:  
  1. Step i-viii of Minimal Gene Specific FLAME  
  8. *splice.signal*  
  9. *shortread*  

- FLAME-GLOW
  ```sh
  ./FLAME.py -I [INPUT.bed] -GTF [Annotation.gtf]
  ```
## Flags
- General FLAME Flags:  
  ```sh
  -I [INPUT.bed]:             #Input file in BED12 format (Required)
  -GTF [Annotation.gtf]:      #Reference annotation file in GTF format (Required)
  --range [int]:              #Variance function window range (Optional, with default = 20)
  -O [string]:                #Output prefix (Optional, with default = "Flame-")
  ```
  
  - Gene Specific FLAME Flags:  
    ```sh
    -G [Gene]:                  #Target gene (Required)
    -R [Reference.fasta]:       #Organism reference assembly sequence in fasta-format for the splice.signal-funciton (Optional)
    -B [Shortread.bam]:         #Short-read RNA sequencing data in bam- or sam-format for the shortread-function (Optional)
    --verbose                   #Optional Flag to output additional files (Optional)
    ```
  
  - FLAME-GLOW Flags:
    ```sh
    --min [int]:                #Minimum read coverage (Optional, with default = 10)
    --ratio [float]:            #Minimum annotation ratio for the FLAME-GLOW module (Optional, with default = 0.25)
    ```
## Test Running FLAME
There is currently only toy data for the ***Gene Specific*** FLAME module.

Both the data as well as instruction for how to run the Gene Specific FLAME module exists under the [`Data`](Data) directory

[//]: <## Output>

[//]: <For more examples and usage, please refer to the [Wiki].>

## Release History
* 0.1.2
    * Updated README.md
    * Added toy data
* 0.1.1
    * FLAME-GLOW mode implemented 
    * Beta-version
* 0.0.1
    * Work in progress
    * Alpha-version

## Meta

Alan Bäckerholm – Alan.baek1@gmail.com

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/marabouboy/](https://github.com/marabouboy/)

## Contributing

1. Fork it (<https://github.com/marabouboy/FLAME/yourfork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
