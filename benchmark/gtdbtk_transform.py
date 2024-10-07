import pandas as pd
import numpy as np
from sys import argv


def importMerged(nodemerged) -> dict:
    merged = pd.read_csv(nodemerged, sep='\t', header=None)
    return merged[[0, 2]].to_dict()


def importNames(names, nodemerged):
    names = pd.read_csv(names, sep='\t', header=None)
    names = names[names[6] == 'scientific name'][[0, 2]]
    names.columns = ['node', 'name']
    node2merged = importMerged(nodemerged)
    for k, v in node2merged.items():
        names.replace(to_replace=k, value=v, inplace=True)
    return  names


if __name__ == '__main__':
    try:
        archaea = pd.read_csv('./gtdbtk.ar122.summary.tsv', sep='\t')
        try:
            bacteria = pd.read_csv('./gtdbtk.bac120.summary.tsv', sep='\t')
            gtdbtk = pd.concat([archaea, bacteria], axis=0)
        except FileNotFoundError:
            gtdbtk = pd.read_csv('./gtdbtk.ar122.summary.tsv', sep='\t')
    except FileNotFoundError:
        gtdbtk = pd.read_csv('./gtdbtk.bac120.summary.tsv', sep='\t')
    # checkm = checkm.merge(checkm_quality[['Bin Id', 'quality',
    #                                       'Completeness',
    #                                       'Contamination']], on='Bin Id')
    # checkm = checkm.replace(np.NaN, 'NOT ASSIGNED')
    # paths = [el.split(';') for el in gtdbtk['classification']]
    lit2taxid = pd.DataFrame(columns=['bins'], data=gtdbtk[['user_genome']])
    lit2taxid['bins'] = gtdbtk['user_genome']
    lit2taxid['paths'] = [el.split(';') for el in gtdbtk['classification']]
    lit2taxid['d__'] = 0
    lit2taxid['p__'] = 0
    lit2taxid['c__'] = 0
    lit2taxid['o__'] = 0
    lit2taxid['f__'] = 0
    lit2taxid['g__'] = 0
    lit2taxid['s__'] = 0
    lit2taxid['class@species'] = False
    lit2taxid['s__name'] = ''
    taxnames = importNames(names='/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/NCBI_taxonomy/from_tar/names.dmp',
                           nodemerged='/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/NCBI_taxonomy/from_tar/merged.dmp')
    taxnames.to_csv('./names.csv')
    ncbi2gtdb = dict()
    i = 0
    for el in ['p__', 'c__', 'o__', 'f__', 'g__', 's__']:
        bact = pd.read_excel('./ncbi_vs_gtdb_r89_bacteria.xlsx', sheet_name=i)
        bact = bact[[bact.columns[0],bact.columns[3]]]
        bact.columns = ['NCBI', 'GTDB']
        arch = pd.read_excel('./ncbi_vs_gtdb_r89_archaea.xlsx', sheet_name=i)
        arch = arch[[arch.columns[0],arch.columns[3]]]
        arch.columns = ['NCBI', 'GTDB']
        ncbi2gtdb[el] = pd.concat([bact, arch], axis=0)
        ncbi2gtdb[el].reset_index(drop=True, inplace=True)
        i += 1
    del bact, arch
    for bin in lit2taxid.bins:
        for name in lit2taxid[lit2taxid.bins == bin]['paths'].iloc[0]:
            n = name[3:]#.replace('_B','').replace('_A', '').replace('_D', '')
            print("printN", n)
            # if n == 'Acinetobacter baumanni':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = taxnames[taxnames['name'] == 'Acinetobacter baumannii'].node.values[0]
            # elif n == 'Bacteroidota':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #     taxnames[taxnames['name'] == 'Bacteroidetes'].node.values[0]
            # elif n == 'Bacteroides vulgatus':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola vulgatus'].node.values[0]
            # elif n == 'Ruminiclostridium thermocellum':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Acetivibrio thermocellus'].node.values[0]
            # elif n == 'Chlorobium phaeoclathratiforme':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Pelodictyon phaeoclathratiforme'].node.values[0]
            # elif n == 'Sulfurihydrogenibium sp000020325':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Pelodictyon phaeoclathratiforme'].node.values[0]
            # elif n == 'Hydrogenobaculum sp000020785':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Hydrogenobaculum sp. Y04AAS1'].node.values[0]
            # elif n == 'Fusobacterium polymorphum':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Fusobacterium nucleatum subsp. polymorphum'].node.values[0]
            # elif n == 'Micromonospora arenicola':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Salinispora arenicola'].node.values[0]
            # elif n == 'Trichormus sp000009705':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Nostoc sp. PCC 7120 = FACHB-418'].node.values[0]
            # elif n == 'Proteiniclasticum sp003514505':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Proteiniclasticum sp.'].node.values[0]
            # elif n == 'Helicobacter_C cinaedi':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Helicobacter_C cinaedi'].node.values[0]
            # elif n == 'Helicobacter pilori_C':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Helicobacter pilori'].node.values[0]
            # elif n == 'Campylobacter showae_C':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Campylobacter showae'].node.values[0]
            # elif n == 'Campylobacter concisus_R':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Campylobacter concisus'].node.values[0]
            # elif n == 'Neisseria_C wadsworthii':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Neisseria wadsworthii'].node.values[0]
            # elif n == 'Neisseria sp000227275' or n == 'Neisseria sp000186165':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Neisseria meningitidis'].node.values[0]
            # elif n == 'Escherichia flexneri':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Escherichia coli'].node.values[0]
            # elif n == 'Bacteroides dorei':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola dorei'].node.values[0]
            # elif n == 'Bacteroides coprophilus':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola coprophilus'].node.values[0]
            # elif n == 'Bacteroides plebeius':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola plebeius'].node.values[0]
            # elif n == 'Bacteroides coprocola':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola coprocola'].node.values[0]
            # elif n == 'Bacteroides sp900066265':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola finegoldii'].node.values[0]
            # elif n == 'Bacteroides sp900066265':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Phocaeicola finegoldii'].node.values[0]
            # elif n == 'Prevotella seregens':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Prevotella dentalis'].node.values[0]
            # elif n == 'Pseudopropionibacterium propionicum':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Arachnia propionica'].node.values[0]
            # elif n == 'Pseudopropionibacterium':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Arachnia'].node.values[0]
            # elif n == 'Propionibacterium humerusii':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Cutibacterium modestum'].node.values[0]
            # elif n == 'Winkia sp002849225':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Winkia neuii'].node.values[0]
            # elif n == 'Pauljensenia sp000185285':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Kytococcus sedentarius'].node.values[0]
            # elif n == 'Rhodococcus hoagii':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Prescottella equi'].node.values[0]
            # elif n == 'Leptotrichia goodfellowii':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Pseudoleptotrichia goodfellowii'].node.values[0]
            # elif n == 'Fusobacterium animalis' or n=='Fusobacterium vincentii':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Fusobacterium nucleatum'].node.values[0]
            # elif n == 'Erysipelatoclostridium spiroforme':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == '[Clostridium] spiroforme'].node.values[0]
            # elif n == 'Absiella dolichum':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Amedibacillus dolichus'].node.values[0]
            # elif n == 'Bacillus thuringiensis_J':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Bacillus thuringiensis'].node.values[0]
            # elif n == 'Lactobacillus ruminis':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Ligilactobacillus ruminis'].node.values[0]
            # elif n == 'Lactobacillus salivarius':
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == 'Ligilactobacillus salivarius'].node.values[0]
            # elif 'Fusobacterium_C' in n:
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == n.replace('_C', '')].node.values[0]
            # elif 'Firmicutes_C' in n:
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == n.replace('_C', '')].node.values[0]
            # elif 'Bacillus_O' in n:
            #     lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #         taxnames[taxnames['name'] == n.replace('_O', '')].node.values[0]
            # elif 'Lactobacillus_H' in n or 'Lactobacillus_G' in n \
            #         or 'Lactobacillus_F' in n or 'Lactobacillus_C' in n:
            #     try:
            #         lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
            #             taxnames[taxnames['name'] == n.replace('_H', '').replace(
            #                 '_G', '').replace('_F', '').replace('_C', '')].node.values[0]
            #     except IndexError:
            #         print(n)
            #
            # else:
            try:
                lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = taxnames[taxnames['name'] == n].node.values[0]
            except IndexError:
                sub = ncbi2gtdb[name[0:3]][ncbi2gtdb[name[0:3]]['GTDB'].apply(lambda x: n in x)]
                print(sub)
                if sub.shape[0] == 0:
                    print('oh no', name)
                    pass
                elif sub.shape[0] == 1:
                    newnm = sub['NCBI'].iloc[0][3:]
                    try:
                        lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
                        taxnames[taxnames['name'] == newnm].node.values[0]
                    except IndexError:
                        print(n, newnm)
                else:
                    print(sub[sub['NCBI'].apply(lambda x: x != name[0:3])]['NCBI'])
                    newnm = sub[sub['NCBI'].apply(lambda x: x != name[0:3])]['NCBI'].iloc[0][3:]
                    try:
                        lit2taxid.loc[lit2taxid.bins == bin, name[0:3]] = \
                        taxnames[taxnames['name'] == newnm].node.values[0]
                    except IndexError:
                        print(n, newnm)


            if name[0:3] == 's__':
                lit2taxid.loc[lit2taxid.bins == bin,'class@species'] = True
                lit2taxid.loc[lit2taxid.bins == bin,'s__name'] = name[3:]




    # newassum = pd.read_csv('/home3/gdefazio/RefSeq_genomes/assembly_summaries_01_10_2020/bacteria_newassum_with_names.csv',
    #                        index_col=0)
    # lvls = ['newtaxid', 'species', 'genus', 'family',
    #         'order', 'class', 'phylum', 'superkingdom']
    # possible_levels = checkm.level.unique()
    # print(possible_levels, lvls.index(possible_levels[0]))
    # paths = {el: newassum[lvls[lvls.index(el):]].drop_duplicates() for el in possible_levels \
    #          if el != 'NOT ASSIGNED'}
    # to_concat = list()
    # for el in possible_levels:
    #     if el != 'NOT ASSIGNED':
    #         to_concat.append(checkm[checkm['level'] == el].merge(paths[el],
    #                                                              left_on='taxid',
    #                                                              right_on=el))
    # to_concat.append(checkm[checkm['level'] == 'NOT ASSIGNED'])
    #
    # total = pd.concat(to_concat, axis=0)
    # if total.shape[0] < checkm.shape[0]:
    #     last = checkm[checkm['Bin Id'].isin(total['Bin Id'])]
    #     paths = newassum[lvls[lvls.index('newtaxid'):]].drop_duplicates()
    #     mgd = last.merge(paths, left_on='taxid', right_on='newtaxid')
    #     total = pd.concat([total, mgd], axis=0)
    #
    # total = total.replace(np.NaN, 0)
    # total = total[['Bin Id', 'taxonomy_lowest_name', 'level', *lvls[1:],
    #                'quality', 'Completeness', 'Contamination']]
    # print('Taxon paths allegation')
    lit2taxid.drop('paths', inplace=True, axis=1)
    lit2taxid.columns = ['bins', 'superkingdom', 'phylum', 'class',
                         'order', 'family', 'genus', 'species',
                         'class@species', 's__name']
    lit2taxid.to_csv('./gtdbtk_assignment_adjusted.csv')


