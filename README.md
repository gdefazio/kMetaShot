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
 conda create --name kmetashot kmetashot -c conda-forge -c gdefazio
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
Before to use kMetaShot you may test the installation typing the following line:

```
kMetaShot_test.py -r /path/to/kMetaShot_reference.h5
```
# USAGE
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
                        to a bin/MAG
  -r , --reference (char)
                        Path to HDF5 kMetaShot reference
  -p , --processes (int)
                        Number of child processes for a Multiprocess parallelism. 
                        Warning: high parallelism <==> high RAM usage
  -o , --out_dir (char)
                        Output directory
  -a , --ass2ref (float)
                        Classification filtering based on ass2ref parameter ranging
                        between 0 and 1. Default 0.
```

kmetashot is also available as Docker container:

```
docker run -it ibiomcnr/kmetashot kMetaShot_classifier_NV.py --help 

                           ################################################
                           #        kMetaShot Classifier Algorithm        #
                           #                  Version 1.0                 #
                           #               Defazio G. et al.              #
                           ################################################

usage: kMetaShot_classifier_NV.py [-h] -b BINS_DIR -r REFERENCE [-a ASS2REF]
                                  -p PROCESSES [-o OUT_DIR]

kMetaShot is able to taxonomically classiy bins/MAGs and long reads by using
an alignment free and k-mer/minimizer based approach.

optional arguments:
  -h, --help            show this help message and exit
  -b BINS_DIR, --bins_dir BINS_DIR
                        Path to directory containing bins or path to multi-
                        fasta file
  -r REFERENCE, --reference REFERENCE
                        Path to HDF5 file containing reference
  -a ASS2REF, --ass2ref ASS2REF
                        Classification filtering based on ass2ref parameter
                        ranging between 0 and 1. Default 0.
  -p PROCESSES, --processes PROCESSES
                        Multiprocess parallelism. Warning: high parallelism
                        <==> high RAM usage
  -o OUT_DIR, --out_dir OUT_DIR
                        Output file path name
```
