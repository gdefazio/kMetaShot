#!/usr/bin/env python3
import os
import h5py
import argparse
import warnings
import argcomplete
import pandas as pd
import multiprocessing as mp
from sys import exit as syexit
from gzip import open as gzopen
from kMetaShot_package.kmer_minimizer_counting_mmh3 import FastaReference, return_time

warnings.filterwarnings('ignore')


def split_options():
    parser = argparse.ArgumentParser(
        description="kMetaShot is able to taxonomically classiy bins/MAGs and long reads by using an alignment free and k-mer/minimizer based approach.",
        prefix_chars="--")
    parser.add_argument("-b", "--bins_dir", type=str,
                        help="Path to directory containing bins or path to multi-fasta file",
                        action="store", required=True)
    parser.add_argument("-r", "--reference", type=str,
                        help="Path to HDF5 file containing reference",
                        action="store", required=True,
                        default=None)
    # parser.add_argument("-y", "--resume_only",
    #                     help="Formulate the resume of bins assignments",
    #                     action="store_true", required=False,
    #                     default=False)
    parser.add_argument("-p", "--processes", type=int,
                        help="Multiprocess parallelism. Warning: high parallelism <==> high RAM usage",
                        action="store", required=True,
                        default=2)
    parser.add_argument("-o", "--out_dir", type=str,
                        help="Output file path name",
                        action="store", required=False,
                        default="./kMetaShot_assignment")
    # parser.add_argument("-c", "--count_ref", type=str,
    #                     help="File with for taxid counts of reference kmers",
    #                     action="store", required=True)
    argcomplete.autocomplete(parser)
    return parser.parse_args()


def adjust_assignmentDF(assignmentDF: pd.DataFrame):
    """
    It prepares the resume table obtained by reference interrogation step
    for the classification step
    :param assignmentDF: raw resume table
    :return: refined resume table
    """
    # i taxid idividuati sono tutti diversi
    # taxid_counts = assignmentDF.taxid.value_counts()
    # se non ci sono taxid rappresentati più di una volta
    # if taxid_counts.shape[0] == assignmentDF.shape[0]:
    md = assignmentDF.merge(assum, on='taxid')
    # print(md)
    # mergio il df delle assegnazioni con assum che
    # presenta i taxid e i relativi path tassonomici
    if md.shape[0] != assignmentDF.shape[0]:
        # se il numero di righe di md è diverso
        # rispetto ad assignment df
        notincluded = assignmentDF[~assignmentDF.taxid.isin(md.taxid)]
        # individuo i taxid che devo includere
        absent_rows = list()
        # includo quel che manca
        for excluded in notincluded.index:
            rowtoadd = assum[assum.genus == notincluded.loc[excluded, 'taxid']]
            rowtoadd = rowtoadd[rowtoadd.index == rowtoadd.index.min()]
            # print(rowtoadd, type(rowtoadd))
            rowtoadd.at[rowtoadd.index.min(), 'taxid'] = 0
            rowtoadd.at[rowtoadd.index.min(), 'species'] = 0
            rowtoadd.at[rowtoadd.index.min(), 'kmer'] = notincluded.loc[excluded, 'kmer']
            absent_rows.append(rowtoadd)
        absent_rows.append(md)
        md = pd.concat(absent_rows, axis=0)
    return md


def assignment_algo(md: pd.DataFrame):
    """
    It is the kMetaShot classification algorithm
    :param md: table resuming number of minimizers of each taxon in q-sequence
    :return: genus, species, strain, ass2ref
    """
    try:
        max_gs = md['genus'].value_counts().index[0]
    except IndexError:
        return 0, 0, 0, 0
    try:
        max_sp = md[(md.genus == max_gs) & (md.species != 0)]['species'].value_counts().index[0]
    except IndexError:
        # print(md['genus'].value_counts())
        # print(md['species'].value_counts())
        return max_gs, 0, 0, 0
    sel_taxids = md[(md.species == max_sp) & (md.taxid != 0)].taxid.unique()
    # ref_kmercount2taxid = pd.read_csv(path_ref_kmercount2taxid, index_col=0)
    ref_taxids = count2taxidref[count2taxidref.taxid.isin(sel_taxids)]
    ass_taxids = pd.DataFrame(md[md.species == max_sp].taxid.value_counts())
    ass_taxids.columns = ['count']
    ass_taxids['taxid'] = ass_taxids.index
    ass_taxids.reset_index(drop=True, inplace=True)
    ass2ref = ass_taxids.merge(ref_taxids, on= 'taxid')
    ass2ref['ratio'] = ass2ref['count_x']/ass2ref['count_y']
    # print(ass_taxids, ref_taxids)
    try:
        strain = ass2ref[ass2ref['ratio'] == ass2ref['ratio'].max()]['taxid'].iloc[0]
        return max_gs, max_sp, strain, ass2ref['ratio'].max()
    except KeyError:
        # print(ass2ref['ratio'].max())
        # print(ass2ref[ass2ref['ratio'] == ass2ref['ratio'].max()].taxid)
        # print(ass2ref)
        return max_gs, max_sp, 0, 0


def reference_importer():
    """
    It is the HDF5 reference File reader.
    :return: kMetaShot reference
    """
    # sostituisci taxid a newtaxid per riadeguare
    ranks = ['newtaxid', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
    assum = pd.read_hdf(reference_path, mode='r', key='new_assemblysummary')
    assum = assum[ranks]
    assum.columns = ['taxid', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
    assum.drop_duplicates(inplace=True)
    return_time('Reference loading ...')
    reference = h5py.File(reference_path, 'r')
    taxid = reference['taxid'][...]
    count2taxidref = pd.read_hdf(reference_path, mode='r', key='count2taxidref',)
    return_time('Reference loading DONE')
    reference.close()
    return taxid, assum, count2taxidref


def single_thread_main(path: str):
    """
    It executes main functions on a single sequence or bin.
    :param path: Path to sequence or bin
    :return: sequence or bin name, genus, species, strain, ass2ref
    """
    return_time("START %s " % path)
    faa = FastaReference(path, cds_ncrna=False)
    faa.seq2bit2kmers(kmer_len, minimizer_len, True)
    bin_out_path = os.path.join(out_dir_deep, '-'.join(path.split('/')[-2:]).split('.f')[0])
    # print(bin_out_path)
    try:
        os.mkdir(bin_out_path)
    except FileExistsError:
        print('%s alreasy exists' % bin_out_path)
    assn = [(k, taxid[k // 8][k % 8]) for k in faa.nphashtable['kmer'] if taxid[k // 8][k % 8] != 0]
    kmer_df = pd.DataFrame(data=assn, columns=['kmer','taxid'])
    adj = adjust_assignmentDF(kmer_df)
    adj.to_csv(os.path.join(bin_out_path, 'kmer2path.csv'),index=False)
    genus, species, strain, ass2ref = assignment_algo(adj)
    return_time('STOP %s' % path)
    return path.split('/')[-1], genus, species, strain, ass2ref


def input_prepare(bins: str) -> list:
    """
    It prepares input in case of FASTA file with multiple sequences.
    It prepares input paths list in case of directory containing bins.
    :param bins: Path to FASTA or path to bins directory
    :return: list of input paths
    """
    todo = list()
    if os.path.isfile(bins):
        # to classify each sequence of a fasta file
        try:
            os.makedirs(os.path.join(out_dir, 'tmp'))
        except FileExistsError:
            pass
        if bins.endswith('.fa') or bins.endswith('.fasta') or bins.endswith('.fna'):
            src = open(bins, 'rt')
        elif bins.endswith('.fa.gz') or bins.endswith('.fasta.gz') or bins.endswith('.fna.gz'):
            src = gzopen(bins, 'rt')
        else:
            raise NotImplementedError('%s is not implemented' % bins)
        c = 0
        for line in src:
            if line.startswith('>'):
                if c == 1:
                    dst.close()
                c = 1
                pathname = os.path.join(out_dir, 'tmp', "%s.fa.gz" %
                                        line[1:-1].replace(
                                            ',', '').replace(
                                            ' ', '_').replace('/', '_'))
                todo.append(pathname)
                dst = gzopen(pathname, 'wt')
                dst.write(line)
            else:
                dst.write(line)
        dst.close()
        src.close()
    else:
        # to classify bins
        for fasta in os.listdir(bins):
            if fasta.endswith('fa') or fasta.endswith('fasta') or fasta.endswith('fna') \
                    or fasta.endswith('fa.gz') or fasta.endswith('fasta.gz') \
                    or fasta.endswith('fna.gz'):
                todo.append(os.path.join(bins, fasta))
            elif os.path.isdir(os.path.join(bins, fasta)):
                for sub in os.listdir(os.path.join(bins, fasta)):
                    pth = os.path.join(bins, fasta, sub)
                    if os.path.isfile(pth):
                        todo.append(pth)
    return todo


if __name__ == '__main__':
    print("")
    print("                           ################################################")
    print("                           #        kMetaShot Classifier Algorithm        #")
    print("                           #                  Version 1.0                 #")
    print("                           #               Defazio G. et al.              #")
    print("                           ################################################")
    print("")
    arguments = split_options()
    return_time('Start Assignment')
    reference_path = arguments.reference
    bins = arguments.bins_dir
    out_dir = arguments.out_dir
    processes = arguments.processes
    kmer_len = 61
    minimizer_len = 31
    print(arguments)
    out_dir_deep = os.path.join(out_dir, 'bins')
    try:
        os.makedirs(out_dir_deep)
    except FileExistsError:
        print("The %s already exists" % out_dir_deep)
    taxid, assum, count2taxidref = reference_importer()
    todo = input_prepare(bins)
    with mp.Pool(processes=processes) as pcs:
        err = pcs.map(single_thread_main,
                      iterable=todo,
                      chunksize=len(todo)//processes)
    pcs.join()
    pcs.close()
    # print(err)
    upper = assum[['genus', 'family', 'order', 'class', 'phylum', 'superkingdom']].drop_duplicates()
    assignment = pd.DataFrame(err, columns=['path', 'genus', 'species', 'strain', 'ass2ref'])
    assignment = assignment[['path', 'ass2ref', 'strain', 'species', 'genus']]
    assignment.columns = ['bin', 'ass2ref', 'taxid', 'species', 'genus']
    zero = assignment[assignment.genus == 0]
    zero['family'] = 0
    zero['order'] = 0
    zero['class'] = 0
    zero['phylum'] = 0
    zero['superkingdom'] = 0
    nonzero = assignment[assignment.genus != 0]
    nonzero = nonzero.merge(upper, on='genus')
    classific = pd.concat([nonzero, zero], axis=0)
    # 1883 specific filter
    classific.loc[(classific.genus == 1883) & (classific.ass2ref < 0.01), [
        'ass2ref', 'taxid', 'species', 'genus','family', 'order', 'class',
        'phylum', 'superkingdom']] = 0
    classific.to_csv(os.path.join(out_dir, 'kMetaShot_classification_resume.csv'))
    return_time('Assignment of bins sequences DONE')
    syexit()
