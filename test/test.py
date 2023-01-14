#!/usr/bin/env python3
import os
import wget
import pandas as pd
import subprocess as sp


if __name__ == '__main__':
    os.mkdir('bins')
    # download genomes to test
    wget.download(url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/921/115/GCF_016921115.1_ASM1692111v1/GCF_016921115.1_ASM1692111v1_genomic.fna.gz',
                  out='./bins/GCF_016921115.1_ASM1692111v1_genomic.fna.gz')
    wget.download(url='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/287/175/GCF_002287175.1_ASM228717v1/GCF_002287175.1_ASM228717v1_genomic.fna.gz',
                  out='./bins/GCF_002287175.1_ASM228717v1_genomic.fna.gz')
    # bins classification test
    sp.run(args=['../kMetaShot_classifier_NV.py',
                 '-b', './bins/',
                 '-r', '../kMetaShot_reference/kMetaShot_bacteria_archaea.h5',
                 '-p', '2',
                 '-o', 'bins_test'], stdout=sp.DEVNULL)

    bins_test = pd.read_csv('./bins_test/kMetaShot_classification_resume.csv', index_col=0)
    if bins_test[bins_test['bin'] == 'GCF_002287175.1_ASM228717v1_genomic.fna.gz'].taxid.iloc[0] == 2161:
        print('Bins test 2161 ... OK')
    if bins_test[bins_test['bin'] == 'GCF_016921115.1_ASM1692111v1_genomic.fna.gz'].taxid.iloc[0] == 1885:
        print('Bins test 1885 ... OK')

    sp.run(args=['../kMetaShot_classifier_NV.py',
                 '-b', './bins/GCF_016921115.1_ASM1692111v1_genomic.fna.gz',
                 '-r', '../kMetaShot_reference/kMetaShot_bacteria_archaea.h5',
                 '-p', '10',
                 '-o', 'contigs_test'], stdout=sp.DEVNULL)
    contigs_test = pd.read_csv('./contigs_test/kMetaShot_classification_resume.csv', index_col=0)
    if contigs_test[contigs_test['bin'] == 'NZ_JAFFZS010000001.1_Streptomyces_actuosus_strain_VRA1_NODE_1_whole_genome_shotgun_sequence.fa.gz'].taxid.iloc[0] == 1885:
        print('Contigs test 1/2 ... OK')
    if contigs_test[contigs_test['bin'] == 'NZ_JAFFZS010000111.1_Streptomyces_actuosus_strain_VRA1_NODE_111_whole_genome_shotgun_sequence.fa.gz'].taxid.iloc[0] == 1885:
        print('Contigs test 2/2 ... OK')
