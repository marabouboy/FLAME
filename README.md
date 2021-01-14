# FLAME: Full Length Adjecency Matrix Enumeration.
FLAME:fire: is a longread splice variant annotation tool. 
It allows for:  
- Filtering and quantification of splice variants and exon connectivity through an [Adjecency Matrix](https://en.wikipedia.org/wiki/Adjacency_matrix)
- Detection of novel exons sites or exon variants through frequency
- Confirmation of said novel exons sites through the use of:
  - Detection of possible adjecent Splice Site Signals
  - Shortread sequencing reads

## Requirements
- Python 3.6
- pysam

## Installation

OS X & Linux:

```
git clone https://github.com/marabouboy/FLAME
```

## Usage example

FLAME runs the Filter, Quantificiation and Novel Exon Detection by default and adheres to the following syntax:
```sh
./FLAME.py -I [INPUT.bed] -GTF [AnnotationReference.gtf] -G [Gene]
```

Flame can be run to confirm the suggested Splice Site Signals and the suggested Exons:
```sh
./FLAME.py -I [INPUT.bed] -GTF [AnnotationReference.gtf] -G [Gene] -R [Reference.fasta] -B [Shortread.bam]/[Shortread.sam]
```

## Flags
```sh
-I [INPUT.bed]:                     #Input file in .bed12 format
-GTF [AnnotationReference.gtf]:     #Reference Annotation file in GTF format
-G [Gene]:                          #Target Gene
--range [int]:                      #Variance Range
-O [string]:                        #Output Prefix
-R [Reference.fasta]:               #Reference in fasta-form to allow for Detection of Adjecent Splice Site Signal 
-B [Shortread.bam]/[Shortread.sam]: #Shortread Sequencing in bam- or sam-format to allow for confirmation of splice site using short read 
```

[//]: <For more examples and usage, please refer to the [Wiki].>

## Development setup

Describe how to install all development dependencies and how to run an automated test-suite of some kind. Potentially do this for multiple platforms.

```
make install
npm test
```

## Release History

* 0.0.1
    * Work in progress

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
