kMetaShot
=========
[![Anaconda-Server Badge](https://anaconda.org/bioconda/kmetashot/badges/downloads.svg)](https://anaconda.org/bioconda/kmetashot)

# Table of content
1. [Introduction](#introduction)  
2. [Install](#install)  
3. [Usage](#usage) 
4. [Citation](#citation)


# INTRODUCTION
The application of 2nd and 3rd generation High Throughput Sequencing (HTS) technologies has deeply reshaped experimental method to investigate microbial communities and obtain a taxonomic and functional profile of the invetigated community. Shotgun Metagenomics allow to quickly obtain a representation of microorganisms genomes characterizing a particular environment. 
In order to obtain a fast e reliable taxonomic classification of microorganisms genomes we present **kMetaShot**, an alignment-free taxonomic classifier based on k-mer/minimizer counting.

# INSTALL
kMetaShot is available through **conda** in _bioconda_ channel. To install it type the following line:  
```
conda create --name kmetashot kmetashot=2.0=pyh7e72e81_1 python=3.9 -c bioconda
```
To activate the environment:
```
conda activate kmetashot
```
## kMetaShot Reference
kMetaShot requires a reference file available at these Zenodo links: 
1. [2nd kMetaShot reference release](https://zenodo.org/records/17375120/files/kMetaShot_bacteria_archaea_2025-05-22.h5) (RefSeq 2025/05/22)
2. [1st kMetaShot reference release](https://zenodo.org/records/17591095/files/kMetaShot_reference.h5) (RefSeq 2022/07/31)

**NEW** kMetaShot reference can be downloaded also from Huggingface in a faster way:
1. [3rd kMetaShot reference release](https://huggingface.co/datasets/gdefazio/kMetaShot_reference_bacteria_archaea/blob/main/compressed_bacteria_archaea_2026-09-10.h5) (RefSeq 2026/09/10)
2. [2nd kMetaShot reference release](https://huggingface.co/datasets/gdefazio/kMetaShot_reference_bacteria_archaea/resolve/main/kMetaShot_bacteria_archaea_2025-05-22.h5) (RefSeq 2025/05/22)
3. [1st kMetaShot reference release](https://huggingface.co/datasets/gdefazio/kMetaShot_reference_bacteria_archaea/resolve/main/kMetaShot_reference.h5) (RefSeq 2022/07/31)

kMetaShot reference represents prokaryotic RefSeq genomes and requires about 22Gb of storage.

## Test
Before to use kMetaShot you may test the installation typing the following line:

```
kMetaShot_test.py -r /path/to/kMetaShot_reference.h5
```

# USAGE
This is the kMetaShot usage.
```
kMetaShot_classifier_NV.py 
                -b bins/
                -r kMetaShot_reference/kMetaShot_bacteria_archaea.h5',
                -p 10
                -o output_dir
                -a 0.1
                
Arguments:
  -h, --help            show this help message and exit
  -b , --bins_dir (char)
                        Path to a directory containing bins fasta files or 
                        path to a multi-fasta file where each header corresponds
                        to a bin/MAG. Files can have .fa, .fasta, .fna, .fa.gz,
                        .fasta.gz, .fna.gz extentions.
  -r , --reference (char)
                        Path to HDF5 kMetaShot reference
  -p , --processes (int)
                        Number of child processes for a Multiprocess parallelism. 
                        Warning: high parallelism <==> high RAM usage
  -o , --out_dir (char)
                        Output directory path
  -a , --ass2ref (float)
                        Classification filtering based on ass2ref parameter ranging
                        between 0 and 1. Default 0. 
                        ass2ref is a ratio between the number of MAG minimizers
                        and the reference minimizers related to the assigned strain
```

kMetaShot is also available as Docker container. It needs `--shm-size=22g` option to properly run as docker container.

```
docker run -it quay.io/biocontainers/kmetashot kMetaShot_classifier_NV.py --help 

```
# kMetaShot as Galaxy tool
kMetaShot has been recently deployed as Galaxy tool in http://usegalaxy.eu server. You can freely 
and easily use it at https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/bgruening/kmetashot/kmetashot/2.0+galaxy2.</br>
</br>
kMetaShot can be also used with a Galaxy instance available at the following</br>
link:</br>
http://212.189.205.125/galaxy/?tool_id=kmetashot&version=latest

# Citation
Giuseppe Defazio, Marco Antonio Tangaro, Graziano Pesole, Bruno Fosso</br>
**kMetaShot: a fast and reliable taxonomy classifier for metagenome-assembled genomes**</br>
Briefings in Bioinformatics, Volume 26, Issue 1, January 2025, bbae680</br>
https://doi.org/10.1093/bib/bbae680
