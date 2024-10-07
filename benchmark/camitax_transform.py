import pandas as pd
import numpy as np
from sys import argv
import warnings

warnings.simplefilter(action='ignore')

if __name__ == '__main__':
    print('Collecting results ./camitax_analysis_resume_paths.csv')
    camitax = pd.read_csv(argv[1], sep='\t')
    camitax.drop('Unnamed: 22', axis=1, inplace=True)
    # print(camitax.columns)
    possible_levels = camitax['taxLvl'].unique()
    # print(possible_levels)
    nodes = pd.read_pickle(
        '/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/NCBI_taxonomy/pandanized/nodes.pkl')
    # camitax = camitax.merge(newassum, left_on='taxID', right_on='taxid', how='left')
    newassum = pd.read_hdf(
        '/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/compressed_clndMashOne4strain_bacteria_archaea.h5',
        key='new_assemblysummary'
    )
    lvls = ['superkingdom', 'phylum', 'class',
            'order', 'family', 'genus', 'species',
            'newtaxid']
    lvls.reverse()
    paths = dict()
    for el in possible_levels:
        if el is not np.NaN:
            paths[el] = newassum[lvls[lvls.index(el):]].drop_duplicates()
    # print(paths.keys())
    to_concat = list()
    for l in possible_levels:
        if l is not np.NaN:
            to_concat.append(camitax[camitax['taxLvl'] == l].merge(paths[l],
                                                                   left_on='taxID',
                                                                   right_on=l))

    total = pd.concat(to_concat, axis=0).replace(np.NaN, 0)
    notincluded = camitax[~camitax.Genome.isin(total.Genome)]

    print(notincluded)
    # notincluded.drop(['superkingdom', 'phylum', 'class',
    #                   'order', 'family', 'genus', 'species'], axis=1, inplace=True)
    for el in notincluded.taxID:
        try:
            rank = nodes[nodes.node == el]['rank'].iloc[0]
        except IndexError:
            print(el)
            break
        if rank in lvls:
            notincluded.loc[notincluded.taxID == el, rank] = el
        parent = nodes[nodes.node == el]['parent'].iloc[0]
        c = 0
        while (rank != 'superkingdom') and (c < 20):
            rank = nodes[nodes.node == parent]['rank'].iloc[0]
            # print(el, rank)
            if rank in lvls:
                notincluded.loc[notincluded.taxID == el, rank] = parent
            parent = nodes[nodes.node == parent]['parent'].iloc[0]
            c += 1

    # if notincluded.shape[0] > 0:
    #     for k in ['species', 'genus', 'family',
    #               'order', 'class', 'phylum','superkingdom']:
    #         notincluded.loc[notincluded.index[0], k] = 0
    total = pd.concat([total, notincluded], axis=0)
    total.to_csv('./camitax_analysis_resume_paths.csv')

