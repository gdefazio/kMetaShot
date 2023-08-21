kMetaShot
=========

# Table of content
1. [Introduction](#introduction)  
2. [Install](#install)  
3. [Usage](#usage)  


# INTRODUCTION
The application of 2nd and 3rd generation High Throughput Sequencing (HTS) technologies has deeply reshaped experimental method to investigate microbial communities and obtain a taxonomic and functional profile of the invetigated community. Shotgun Metagenomics allow to quickly obtain a representation of microorganisms genomes characterizing a particular environment. 
In order to obtain a fast e reliable taxonomic classification of microorganisms genomes we present **kMetaShot**, an alignment-free taxonomic classifier based on k-mer/minimizer counting.

# INSTALL
kMetaShot is available through **conda**. To install it type the following line:  
```
 conda create —name kmetashot kmetashot -c conda-forge -c mtangaro
```
To activate the environment:
```
conda activate kmetashot
```

kMetaShot requires a reference file available at this [link](http://srv00.recas.ba.infn.it/webshare/brunofosso/kMetaShot_reference.h5). It requires about 22Gb of storage.  
To download it you can simply use `wget`:
```
wget http://srv00.recas.ba.infn.it/webshare/brunofosso/kMetaShot_reference.h5
```

# USAGE
```
kMetaShot_classifier_NV.py 
                -b bins/
                -r kMetaShot_reference/kMetaShot_bacteria_archaea.h5',
                -p 10
                -o output_dir
                
Arguments:
  -h, --help            show this help message and exit
  -b , --bins_dir (char)
                        Path to a directory containing bins fasta files or 
                        path to a multi-fasta file where each header corresponds
                        to a bin/MAG
  -r , --reference (char)
                        Path to HDF5 kMetaShot reference
  -p , --processes (int)
                        Number of child processes for a Multiprocess parallelism. 
                        Warning: high parallelism <==> high RAM usage
  -o , --out_dir (char)
                        Output directory
```
