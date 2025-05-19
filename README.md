kMetaShot
=========

# Table of content
1. [Introduction](#introduction)  
2. [Install](#install)  
3. [Usage](#usage) 
4. [Citation](#citation)


# INTRODUCTION
The application of 2nd and 3rd generation High Throughput Sequencing (HTS) technologies has deeply reshaped experimental method to investigate microbial communities and obtain a taxonomic and functional profile of the invetigated community. Shotgun Metagenomics allow to quickly obtain a representation of microorganisms genomes characterizing a particular environment. 
In order to obtain a fast e reliable taxonomic classification of microorganisms genomes we present **kMetaShot**, an alignment-free taxonomic classifier based on k-mer/minimizer counting.

# INSTALL
kMetaShot is available through **conda**. To install it type the following line:  
```
 conda create --name kmetashot kmetashot=2.0=py_0 -c conda-forge -c gdefazio
```
To activate the environment:
```
conda activate kmetashot
```

kMetaShot requires a reference file available at this [link](http://srv00.recas.ba.infn.it/webshare/brunofosso/kMetaShot_reference.h5). It requires about 22Gb of storage.
It can be used also the following command:
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

kmetashot is also available as Docker container:

```
docker run -it ibiomcnr/kmetashot kMetaShot_classifier_NV.py --help 

```
# Galaxy
kMetaShot can be also used by employing a Galaxy instance available at the following</br>
link:</br>
http://212.189.202.40.cloud.ba.infn.it/galaxy/?tool_id=kmetashot_pipeline&version=latest

# Citation
Giuseppe Defazio, Marco Antonio Tangaro, Graziano Pesole, Bruno Fosso</br>
**kMetaShot: a fast and reliable taxonomy classifier for metagenome-assembled genomes**</br>
Briefings in Bioinformatics, Volume 26, Issue 1, January 2025, bbae680</br>
https://doi.org/10.1093/bib/bbae680
