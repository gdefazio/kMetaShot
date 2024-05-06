import pandas as pd
from sys import argv
from os.path import join as pjoin
import warnings

warnings.simplefilter(action='ignore')

sample_dir = argv[1]
metadata_pth = argv[2]
sample_name = "_".join(sample_dir.split('/')[-2].split('_')[-2:])
sample_mapping = pd.read_csv(pjoin(sample_dir, 'gsa_mapping.tsv.gz'),
                             sep='\t')
genomeid = list(sample_mapping['genome_id'].unique())
metadata = pd.read_csv(metadata_pth, sep='\t')
sample_comp = metadata[metadata['genome_ID'].isin(genomeid)]
# sample_comp.to_csv('./sample12_composition.csv')

assum = pd.read_hdf('/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/clndMashOne4strain_bacteria_archaea.h5','new_assemblysummary')

taxids = ['superkingdom', 'phylum', 'class', 'order',
          'family', 'genus', 'species', 'taxid']
taxid = assum[taxids].drop_duplicates()
species = assum[taxids[:-1]].drop_duplicates()
# genus = assum[taxids[:-2]].drop_duplicates()

# a = pd.read_csv('./sample12_composition.csv', index_col=0)

comp_tax_lvl = sample_comp.merge(taxid, left_on='NCBI_ID',right_on='taxid')
comp_spec_lvl = sample_comp.merge(species, left_on='NCBI_ID',right_on='species')

inref = pd.concat([comp_tax_lvl,
                   comp_spec_lvl[~comp_spec_lvl.NCBI_ID.isin(comp_tax_lvl.NCBI_ID)]],
                  axis=0)

inref['inref'] = True
taxlin = pd.read_pickle('/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/NCBI_taxonomy/pandanized/taxidlineage.pkl')

outref = sample_comp[~sample_comp.NCBI_ID.isin(inref.NCBI_ID)]

nodes = pd.read_pickle('/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/NCBI_taxonomy/pandanized/nodes.pkl')


lvls = ['superkingdom', 'phylum', 'class',
            'order', 'family', 'genus', 'species',
            'taxid']
for el in outref.NCBI_ID:
    if el == 1834200:
        el = 1796646
    elif el == 46170:
        el = 1280
    try:
        rank = nodes[nodes.node == el]['rank'].iloc[0]
        if rank in lvls:
            outref.loc[outref.NCBI_ID == el, rank] = el
        parent = nodes[nodes.node == el]['parent'].iloc[0]
        c = 0
        while (rank != 'superkingdom') and (c < 20):
            rank = nodes[nodes.node == parent]['rank'].iloc[0]
            # print(el, rank)
            if rank in lvls:
                outref.loc[outref.NCBI_ID == el, rank] = parent
            parent = nodes[nodes.node == parent]['parent'].iloc[0]
            c += 1
    except Exception:
        print(el)


outref['inref'] = False
composition = pd.concat([inref, outref], axis=0)
composition.to_csv(pjoin(sample_dir, '%s_composition_fullpath.csv' % sample_name))
print(sample_name)