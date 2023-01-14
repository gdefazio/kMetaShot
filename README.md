# kMetaShot

The application of 2nd and 3rd generation High Throughput Sequencing (HTS) technologies has deeply reshaped experimental method to investigate microbial communities to obtain a compositional landscape. Shotgun Metagenomics allow to quickly obtain a representation of microorganisms genomes characterizing a particular environment. Then, in order to obtain a fast e reliable taxonomic classification of microorganisms genomes we present kMetaShot, an alignment-free taxonomic classifier based on k-mer/minimizer counting.

### USAGE
`{kMetaShot_classifier_NV.py 
                -b bins/
                -r kMetaShot_reference/kMetaShot_bacteria_archaea.h5',
                -p 10
                -o output_dir
}`