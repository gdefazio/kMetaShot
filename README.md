# kMetaShot

The application of 2nd and 3rd generation High Throughput Sequencing (HTS) technologies has deeply reshaped experimental method to investigate microbial communities and obtain a compositional landscape. Shotgun Metagenomics allow to quickly obtain a representation of microorganisms genomes characterizing a particular environment. Then, in order to obtain a fast e reliable taxonomic classification of microorganisms genomes we present kMetaShot, an alignment-free taxonomic classifier based on k-mer/minimizer counting.

### USAGE
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