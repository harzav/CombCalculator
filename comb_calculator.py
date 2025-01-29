#!/usr/bin/env python3

"""
example of calling the script:
python3.10 /dir/CombCalculator/comb_calculator.py \
--input /dir/CombCalculator/example_inputs/uniprot_mapped_disorder_TEST.csv \
--output /dir/callable_output_TEST/ \
--max_combinations 1000000 \
--threads 4 \
--batch_size 5 \
--additional /dir/CombCalculator/example_inputs/combinationsTEST.csv

"""
import argparse
import itertools
import os
from pathlib import Path
import time
import functools
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), 'pygosemsim'))
from pygosemsim import term_set
from pygosemsim import annotation
from pygosemsim import similarity
from pygosemsim import graph
from goatools import obo_parser
import wget
from Bio import Align
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import joblib


BASE_DIR = Path(__file__).resolve().parent
DATASETS_DIR = BASE_DIR / "Datasets" 

# Ensure directories exist
DATASETS_DIR.mkdir(parents=True, exist_ok=True)

#-----------------------------------------------------------------------------------------------------------------------#
# Load already created datasets for checking existence of ppis
#-----------------------------------------------------------------------------------------------------------------------#

def normalize_pairs(df):
    """Normalize the uidA and uidB columns for consistency."""
    df[['uidA', 'uidB']] = pd.DataFrame(df[['uidA', 'uidB']].apply(sorted, axis=1).to_list(), index=df.index)
    return df

def read_and_normalize_input(input_path):
    """Read the input dataset and normalize the uid columns."""
    df = pd.read_csv(input_path)
    #df = normalize_pairs(df)
    return df

def create_df_to_be_excluded(additional_df_path=None):
    """Create a DataFrame to be excluded from combinations."""
    # Load and normalize data as previously implemented
    dfTR = pd.read_csv(DATASETS_DIR / "taste_receptor_comb.csv", usecols=['uidA', 'uidB'])
    dfTRAIN = pd.read_csv(DATASETS_DIR / "FINAL_DATASET.csv", usecols=['uidA', 'uidB'])
    dfTEST = pd.read_csv(DATASETS_DIR / "FINAL_TEST_DS.csv", usecols=['uidA', 'uidB'])
    df_db = pd.read_csv(DATASETS_DIR / "ppis_included.csv", usecols=['uidA', 'uidB'])

    # Normalize the datasets
    dfTR = normalize_pairs(dfTR)
    dfTRAIN = normalize_pairs(dfTRAIN)
    dfTEST = normalize_pairs(dfTEST)
    df_db = normalize_pairs(df_db)

    # Combine datasets
    combined_df = pd.concat([dfTR, dfTRAIN, dfTEST, df_db], ignore_index=True)

    if additional_df_path:
        user_df = pd.read_csv(additional_df_path)
        if 'uidA' not in user_df or 'uidB' not in user_df:
            raise ValueError(f"Dataset at {additional_df_path} must contain 'uidA' and 'uidB' columns.")
        user_df = normalize_pairs(user_df)
        combined_df = pd.concat([combined_df, user_df], ignore_index=True)

    # Drop duplicates and return
    df_to_be_excluded = combined_df.drop_duplicates(subset=['uidA', 'uidB'], keep='first')
    return df_to_be_excluded

#-----------------------------------------------------------------------------------------------------------------------#
# Load datasets for feature calculation
#-----------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------HOMOLOGY DATASETS-------------------------------------------------------------#

homologene_data=pd.read_csv(DATASETS_DIR / "homology files" / "homologene_data.csv")
homologene_data['protein_accession'] = homologene_data['protein_accession'].str.split('.').str[0]

#import mouse
map_mouse=pd.read_csv(DATASETS_DIR / "homology files" / "map_mouse.csv")
map_mouse['Protein_Accession'] = map_mouse['Protein_Accession'].str.split('.').str[0]
mouse_ppi=pd.read_csv(DATASETS_DIR / "homology files" / "mouse_ppi.csv")
mouse_ppi.rename(columns={'ID interactor A':'uniprotid_A', 'ID interactor B':'uniprotid_B'}, inplace=True)

#import drosophila
map_drosophila=pd.read_csv(DATASETS_DIR / "homology files" / "map_drosophila.csv")
map_drosophila['Protein_Accession'] = map_drosophila['Protein_Accession'].str.split('.').str[0]
dros_ppi=pd.read_csv(DATASETS_DIR / "homology files" / "drosophila_ppi.csv")
dros_ppi.rename(columns={'ID interactor A':'uniprotid_A','ID interactor B':'uniprotid_B'}, inplace=True)
#import yeast
map_yeast=pd.read_csv(DATASETS_DIR / "homology files" / "map_yeast.csv")
map_yeast['Protein_Accession'] = map_yeast['Protein_Accession'].str.split('.').str[0]
yeast_ppi=pd.read_csv(DATASETS_DIR / "homology files" / "yeast_ppi.csv")
yeast_ppi.rename(columns={'ID interactor A':'uniprotid_A', 'ID interactor B':'uniprotid_B'}, inplace=True)
#import Ecoli
map_ecoli=pd.read_csv(DATASETS_DIR / "homology files" / "map_Ecoli.csv")
map_ecoli['Protein_Accession'] = map_ecoli['Protein_Accession'].str.split('.').str[0]
ecoli_ppi=pd.read_csv(DATASETS_DIR / "homology files" / "Ecoli_ppi.csv")
ecoli_ppi.rename(columns={'ID interactor A':'uniprotid_A', 'ID interactor B':'uniprotid_B'}, inplace=True)
mint_ds = pd.read_csv(DATASETS_DIR / "homology files" / "human_MINT_curated.csv")
dip_ds = pd.read_csv(DATASETS_DIR / "homology files" / "human_DIP_curated.csv")
apid_ds = pd.read_csv(DATASETS_DIR / "homology files" / "human_APID_curated.csv")
biogrid_ds = pd.read_csv(DATASETS_DIR / "homology files" / "human_BIOGRID_curated.csv")

#-----------------------------------------SL DATASETS-------------------------------------------------------------#

map_ds = pd.read_csv(DATASETS_DIR / "SL files" / "map_uniprot_human.csv", usecols=['ID', 'Protein_Accession'])
es_ds = pd.read_table(DATASETS_DIR / "SL files" / "eSLDB_Homo_sapiens.txt", delimiter='\t', usecols=['Experimental annotation', 'SwissProt entry'])
es_ds = es_ds[~es_ds["Experimental annotation"].isnull()]
es_ds = es_ds[~es_ds["SwissProt entry"].isnull()]
es_ds.rename(columns={'Experimental annotation': 'Experimental_annotation', 'SwissProt entry': 'Protein_Accession'}, inplace=True)
es_ds['Experimental_annotation'] = es_ds['Experimental_annotation'].str.split(',')

#-----------------------------------------GD DATASETS-------------------------------------------------------------#

map_test=pd.read_csv(DATASETS_DIR / 'spearman files' / 'map.csv')
map_test.rename(columns={'UniProtKB-AC':'UniprotID'}, inplace=True)
df1 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS181.soft', delimiter="\t", skiprows=[i for i in range(306,)])
df2 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS531.soft', delimiter="\t", skiprows=[i for i in range(210,)])
df3 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS807.soft', delimiter="\t", skiprows=[i for i in range(98,)], comment='!')
df4 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS841.soft', delimiter="\t", skiprows=[i for i in range(85,)], comment='!')
df5 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS806.soft', delimiter="\t", skiprows=[i for i in range(97,)], comment='!')
df6 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS987.soft', delimiter="\t", skiprows=[i for i in range(78,)])
df7 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS2855.soft', delimiter="\t", skiprows=[i for i in range(202,)])
df8 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS1088.soft', delimiter="\t", skiprows=[i for i in range(58,)])
df9 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS843.soft', delimiter="\t", skiprows=[i for i in range(113,)], comment='!')
df10 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS1402.soft', delimiter="\t", skiprows=[i for i in range(262,)])
df11 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS1085.soft', delimiter="\t", skiprows=[i for i in range(162,)], comment='!')
df12 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS534.soft', delimiter="\t", skiprows=[i for i in range(117,)])
df13 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS3257.soft', delimiter="\t", skiprows=[i for i in range(199,)])
df14 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS651.soft', delimiter="\t", skiprows=[i for i in range(78,)])
df15 = pd.read_table(DATASETS_DIR / 'spearman files' / 'GDS596.soft', delimiter="\t", skiprows=[i for i in range(580,)])
GDtotal=[df1, df2, df3 , df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15]
GD_test=[]
for i in GDtotal:
    i_numeric = i.iloc[2:].select_dtypes(include=[np.number])
    i = i.fillna(i_numeric.mean())
    i=i.drop(columns=['ID_REF'])
    i=i.groupby('IDENTIFIER').mean().reset_index()
    GD_test.append(i)


#-----------------------------------------RNASEQ EXPRESSION DATASETS-------------------------------------------------------------#

expression_ds=pd.read_csv(DATASETS_DIR / "extra rna seq columns" / "GSE227375_norm_counts_STAR.csv")
expression_ds.columns = expression_ds.columns.str.replace('Unnamed: 0','ENS_ID')
map_ds=pd.read_csv(DATASETS_DIR / "extra rna seq columns" / "UID_to_ENSBL_mapping.csv", usecols=['uid', 'ID'])
map_ds['ID'] = map_ds['ID'].str.split('.').str[0]
expression_ds_2=pd.read_csv(DATASETS_DIR / "extra rna seq columns" / "GSE228702_adjusted_expression_Cellcounts_granulatorAbis0_nnls.csv")
expression_ds_2.columns = expression_ds_2.columns.str.replace('Unnamed: 0','ENS_ID')
map_ds_2=pd.read_csv(DATASETS_DIR / "extra rna seq columns" / "UID_to_GeneCards_mapping.csv", usecols=['uid', 'ID'])
map_ds_2['ID'] = map_ds_2['ID'].str.split('.').str[0]

#-----------------------------------------------------------------------------------------------------------------------#
# Feature calculation functions
#-----------------------------------------------------------------------------------------------------------------------#

def interpro_calculation (uid_a, uid_b):
    def fetch_pfam_ids(prot):
        pfam_ids = {}
        url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/UniProt/{prot}/"
        try:
            response = requests.get(url)
            if response.status_code == 200:
                payload = response.json()
                if 'results' in payload and isinstance(payload['results'], list):
                    ipr_list = [item['metadata']['accession'] for item in payload['results'] if item.get('metadata', {}).get('source_database') == 'pfam']
                    pfam_ids[prot] = ipr_list if ipr_list else 'NaN'
                else:
                    pfam_ids[prot] = 'NaN'
            else:
                pfam_ids[prot] = 'NaN'
        except (requests.RequestException, ValueError):
            pfam_ids[prot] = 'NaN'
        return pfam_ids


    pfam_A = fetch_pfam_ids(uid_a)
    pfam_B = fetch_pfam_ids(uid_b)

    pfam_pair_list = [pfam_A.get(uid_a, 'NaN'), pfam_B.get(uid_b, 'NaN')]


    with open(DATASETS_DIR / "3did_flat_Mar_4_2021.dat") as file:
        interact_list = [x.split("\t")[3:5] for x in file.readlines() if x.startswith('#=ID')]
        interact_list = [[p_A.split(".")[0], p_B.split(".")[0]] for p_A, p_B in interact_list]
        for sublist in interact_list:
            sublist[0] = sublist[0].lstrip(' (')
    


    result_list = []
    for interaction in interact_list:
        if (interaction[0] in pfam_pair_list[0] and interaction[1] in pfam_pair_list[1]) or (interaction[0] in pfam_pair_list[1] and interaction[1] in pfam_pair_list[0]):
            result_list.append(1)
        elif 'NaN' in pfam_pair_list:
            result_list.append('NaN')
        else:
            result_list.append(0)

    split_sp_lists = [result_list[x:x+len(interact_list)] for x in range(0, len(result_list), len(interact_list))]
    final_list = [1 if 1 in sublist else ('NaN' if 'NaN' in sublist else 0) for sublist in split_sp_lists]

    return final_list[0]

def e_value_calc (seq_a, seq_b):
    aligner = Align.PairwiseAligner(match_score=1.0)
    if (seq_a != 'NaN') and (seq_b != 'NaN'):
        score = aligner.score(seq_a, seq_b)
        mean_len = (len(seq_a) + len(seq_b)) / 2
        e_val = mean_len * (2.7**(-score))
        return e_val
    else:
        return 'NaN'

def extra_feature_calc (seq_a, seq_b):
    
    sequences1 = [[seq_a, seq_b]]
    sequences = [[i[0].replace('U', ''), i[1].replace('U', '')] for i in sequences1]
    sequences = [[i[0].replace('X', ''), i[1].replace('X', '')] for i in sequences]
    sequences = [[i[0].replace(' ', ''), i[1].replace(' ', '')] for i in sequences]
    sequences = [[i[0].replace('\n', ''), i[1].replace('\n', '')] for i in sequences]


    def get_features(X):
        amino_acids = ['A', 'L', 'I', 'M', 'F', 'V', 'S', 'P', 'T', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'G']
        aaper = X.get_amino_acids_percent()
        features = [aaper[aa] for aa in amino_acids]

        mw = X.molecular_weight()
        arom = X.aromaticity()
        insta = X.instability_index()
        fraction = X.secondary_structure_fraction()
        helix_fraction, turn_fraction, sheet_fraction = fraction[0], fraction[1], fraction[2]
        extinct = X.molar_extinction_coefficient()
        cys_reduced, cys_residues = extinct[0], extinct[1]
        gravy = X.gravy(scale='KyteDoolitle')
        ph_charge = X.charge_at_pH(7)

        features.extend([mw, arom, insta, helix_fraction, turn_fraction, sheet_fraction,
                        cys_reduced, cys_residues, gravy, ph_charge])
        return features

    # Calculate features for each sequence pair
    feature_list = []
    for seq in sequences:
        X0 = ProteinAnalysis(seq[0])
        seq0 = get_features(X0)
        X1 = ProteinAnalysis(seq[1])
        seq1 = get_features(X1)
        substr_list = [abs(x - y) for x, y in zip(seq0, seq1)]
        abs_list = [abs(val) for val in substr_list]
        feature_list.append(abs_list)

    return feature_list[0]

def go_similarity_code(uid_a, uid_b):
    
    G = graph.from_resource("go-basic")
    similarity.precalc_lower_bounds(G)
    annot = annotation.from_resource("goa_human")

    test_list = [[uid_a, uid_b]]

    annot_list_full = []
    for list1 in test_list:
        annot_list = []
        for i in list1:
            try:
                annot1 = str(annot[i]["annotation"].keys())
                annot2 = annot1.strip('dict_keys([])')
                annot2 = annot2.strip("'")
                annot3 = list(annot2.split("', '"))
                annot_list.append(annot3)
            except KeyError:
                annot_list.append("NaN")
        annot_list_full.append(annot_list)

    go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
    data_folder = DATASETS_DIR / 'data'  # Save the data to this directory

    # Check if the path exists as a file
    if os.path.isfile(data_folder):
        raise Exception('Data path (' + data_folder + ') exists as a file. '
                        'Please rename, remove, or change the desired location of the data path.')

    # Check if the directory already exists
    if not os.path.isdir(data_folder):
        try:
            os.mkdir(data_folder)
        except OSError as e:
            if e.errno != 17:  # Ignore "directory exists" error
                raise e

    # Download the GO OBO file if it doesn't exist
    go_obo_file = os.path.join(data_folder, 'go-basic.obo')
    if not os.path.isfile(go_obo_file):
        go_obo = wget.download(go_obo_url, go_obo_file)
    else:
        go_obo = go_obo_file

    # Load the GO DAG
    go = obo_parser.GODag(go_obo)

    def add_similarity_to_test(similarity_name, similarity_df_header):
        go_sim_list=[]
        for annot_list in annot_list_full:
            def filter_func(filter_term):
                BP_list_A=[]
                BP_list_B=[]
                BP_list=[BP_list_A,BP_list_B]
                if annot_list[0]!='NaN':
                    for i in annot_list[0]:
                        go_id = i
                        try:
                            go_term = go[go_id]
                            if go_term.namespace == filter_term:
                                BP_list_A.append(i)
                        except KeyError:
                            continue
                else:
                    BP_list_A.append('NaN')
                if annot_list[1]!='NaN':
                    for i in annot_list[1]:
                        go_id = i
                        try:
                            go_term = go[go_id]
                            if go_term.namespace == filter_term:
                                BP_list_B.append(i)
                        except KeyError:
                            continue
                else:
                    BP_list_B.append('NaN')
                if ((BP_list_A != 'NaN') & (BP_list_A != 'NaN')):
                    sf = functools.partial(term_set.sim_func, G, similarity.lin)
                    go_sim=term_set.sim_bma(BP_list[0], BP_list[1], sf)
                    return go_sim
                else:
                    return 'NaN'

            bp_sim= filter_func(similarity_name)
            go_sim_list.append(bp_sim)

        return go_sim_list

    bp_list= add_similarity_to_test('biological_process', 'BP_similarity')
    mf_list= add_similarity_to_test('molecular_function', 'MF_similarity')
    cc_list=add_similarity_to_test('cellular_component', 'CC_similarity')

    return bp_list, mf_list, cc_list



def homology_code (refseq_id_a, refseq_id_b, uid_a, uid_b):

    test_list =[[refseq_id_a, refseq_id_b]]

    homol_list=[]
    for pair in test_list:
        hom1=homologene_data[homologene_data["protein_accession"].isin(pair) == True]
        homol_list.append(hom1)
    homol_list2=[]
    for i in homol_list:
        hom1=homologene_data[homologene_data["HID"].isin(i['HID']) == True]
        homol_list2.append(hom1)

    def add_column_function(map_ds, ppi_ds,taxid_no):
        homol_list_filtered=[]
        for i in homol_list2:
            hom1=i.loc[i['taxid']==taxid_no]
            homol_list_filtered.append(hom1)

        merge_list=[]
        for ds in homol_list_filtered:
            mrg= pd.merge ( ds,map_ds, left_on='protein_accession', right_on='Protein_Accession')
            merge_list.append(mrg)

        pair_list=[]
        nan_list=['NaN','NaN']
        for i in merge_list:
            if i['HID'].nunique()==2:
                sublist=list(i['UniprotID'])
                pair_list.append(sublist)
            else:
                pair_list.append(nan_list)
        ppi_list=[]
        for index, rows in ppi_ds.iterrows():
            my_list =[rows.uniprotid_A, rows.uniprotid_B]
            ppi_list.append(my_list)
        
        value_list=[]
        for i in pair_list:
            for j in ppi_list:
                if (len(i)>=2)& (set(i)==set(j)):
                    value_list.append(1)
                elif 'NaN' in i:
                    value_list.append('NaN') 
                else:
                    value_list.append(0)

        split_sp_lists = [value_list[x:x+len(ppi_list)] for x in range(0, len(value_list), len(ppi_list))]
        final_list=[]
        for i in split_sp_lists:
            if 1 in i:
                final_list.append(1)
            elif 'NaN' in i:
                final_list.append('NaN')
            else:
                final_list.append(0)

        return final_list[0]

    mouse_list=add_column_function(map_mouse, mouse_ppi, 10090)
    dros_list=add_column_function(map_drosophila, dros_ppi, 7227)
    yeast_list=add_column_function(map_yeast, yeast_ppi, 559292)
    ecoli_list= add_column_function(map_ecoli, ecoli_ppi, 83333)

    def database_search(database_dataset):
        pair_list=[[uid_a, uid_b]]

        ppi_list=[]
        for index, rows in database_dataset.iterrows():
            my_list =[rows.uidA, rows.uidB]
            ppi_list.append(my_list)

        value_list=[]
        for i in pair_list:
            for j in ppi_list:
                if set(i)==set(j):
                    value_list.append(1)
                else:
                    value_list.append(0)

        split_sp_lists = [value_list[x:x+len(ppi_list)] for x in range(0, len(value_list), len(ppi_list))]
        final_list=[]
        for i in split_sp_lists:
            if 1 in i:
                final_list.append(1)
            else:
                final_list.append(0)

        return final_list[0]

    mint=database_search(mint_ds)
    dip=database_search(dip_ds)
    apid=database_search(apid_ds)
    biogrid=database_search( biogrid_ds)

    return mouse_list, dros_list, yeast_list, ecoli_list, mint, dip, apid, biogrid



def sl_calculation (uid_a, uid_b):

    pair_set = {(uid_a, uid_b)}

    values_list = []
    for pair in pair_set:
        if pair[0] == pair[1]:
            values_list.append(1)
        else:
            match = map_ds[map_ds["ID"].isin(pair)]
            if len(match) > 0:
                merged = pd.merge(match, es_ds, on='Protein_Accession')
                annotations = merged['Experimental_annotation'].tolist()
                values_list.append(annotations if len(annotations) > 1 else 'NaN')
            else:
                values_list.append('NaN')

    def compare_annotations(annotation_list):
        if annotation_list == 1:
            return 1
        elif annotation_list == 'NaN':
            return 'NaN'
        else:
            annotations_A = set([annotation.strip() for annotation in annotation_list[0]])
            annotations_B = set([annotation.strip() for annotation in annotation_list[1]])
            if annotations_A.intersection(annotations_B):
                return 1
            else:
                return 0

    def process_annotations(val):
        if val == 1:
            return 1
        elif val == 'NaN':
            return 'NaN'
        else:
            return compare_annotations(val)

    final_list = [process_annotations(val) for val in values_list]
    return final_list[0]





def gd_datasets_calculation (uid_a, uid_b):
        
    uniprot_list =[[uid_a, uid_b]]

    map_list=[]
    for row in uniprot_list:
        mapds1=map_test[map_test["UniprotID"].isin(row) == True]
        map_list.append(mapds1)
    
    pair_list=[]
    for i in GD_test:
        for ds in map_list:
            GD1=i[i["IDENTIFIER"].isin(ds["ID"]) == True]
            pair_list.append(GD1)

    sp_list=[]
    for pair in pair_list:
        if len(pair.index)==2:
            pair=pair.drop(columns=['IDENTIFIER'])
            pair2=pair.transpose()
            rho, p= spearmanr(pair2)
            sp_list.append(rho)
        else:
            sp_list.append('nan')
    
    return sp_list


def rna_dataset_1_calculation (uid_a, uid_b):
    data = {'uidA': [uid_a],
        'uidB': [uid_b]}
    test_uids= pd.DataFrame(data)

    test_mapped= pd.merge(test_uids, map_ds, left_on='uidA', right_on='uid', how='left')
    test_mapped.rename(columns={'ID':'ID_A'}, inplace=True)
    test_mapped= pd.merge(test_mapped, map_ds, left_on='uidB', right_on='uid', how='left')
    test_mapped.rename(columns={'ID':'ID_B'}, inplace=True)

    #---------------- A-------------------------------------#
    map_merge_A= pd.merge(test_mapped, expression_ds, left_on='ID_A', right_on='ENS_ID', how='left')
    del map_merge_A['uid_x']
    del map_merge_A['uid_y']
    #------ calculate average between expressions lists of different RNAs of the same protein
    group_map_merge_A=map_merge_A.groupby(['ENS_ID'], sort=False).mean(numeric_only=True)
    group_map_merge_A['expression_list_A']= group_map_merge_A.loc[:, group_map_merge_A.columns != 'ENS_ID'].values.tolist()

    map_merge_A=map_merge_A[['uidA', 'uidB','ENS_ID']]
    merge_A=pd.merge(map_merge_A, group_map_merge_A,on='ENS_ID', how='left')
    merge_A=merge_A[['uidA', 'uidB','expression_list_A']]


    #---------------------- B--------------------------#
    map_merge_B= pd.merge(test_mapped, expression_ds, left_on='ID_B', right_on='ENS_ID', how='left')
    del map_merge_B['uid_x']
    del map_merge_B['uid_y']
    group_map_merge_B=map_merge_B.groupby(['ENS_ID'], sort=False).mean(numeric_only=True)
    group_map_merge_B['expression_list_B']= group_map_merge_B.loc[:, group_map_merge_B.columns != 'ENS_ID'].values.tolist()

    map_merge_B=map_merge_B[['uidA', 'uidB','ENS_ID']]
    merge_B=pd.merge(map_merge_B, group_map_merge_B,on='ENS_ID', how='left')
    merge_B=merge_B[['uidA', 'uidB','expression_list_B']]


    # DROP DUPLICATES AND! NANs
    #----------------------------------------------------------#
    merge_A = merge_A.drop_duplicates(subset=['uidA','uidB'])
    merge_B = merge_B.drop_duplicates(subset=['uidA','uidB'])
    merge_B =merge_B.drop(['uidA', 'uidB'], axis=1)
    #----------------------------------------------------------#

    concatds=pd.concat([merge_A, merge_B], axis=1)

    exp_list =[]
    for index, rows in concatds.iterrows():
        my_list =[rows.expression_list_A, rows.expression_list_B]
        exp_list.append(my_list)

    calc_list=[]
    for pair in exp_list:
        if (str(pair[0])=='nan')|(str(pair[1])=='nan'):
            calc_list.append('NaN')
        else:
            
            rho, p= spearmanr(pair[0], pair[1])
            calc_list.append(rho)

    return calc_list[0]



def rna_dataset_2_calculation (uid_a, uid_b):
    data = {'uidA': [uid_a],
        'uidB': [uid_b]}
    test_uids= pd.DataFrame(data)

    test_mapped= pd.merge(test_uids, map_ds_2, left_on='uidA', right_on='uid', how='left')
    test_mapped.rename(columns={'ID':'ID_A'}, inplace=True)
    test_mapped= pd.merge(test_mapped, map_ds_2, left_on='uidB', right_on='uid', how='left')
    test_mapped.rename(columns={'ID':'ID_B'}, inplace=True)

    #---------------- A-------------------------------------#
    map_merge_A= pd.merge(test_mapped, expression_ds_2, left_on='ID_A', right_on='ENS_ID', how='left')
    del map_merge_A['uid_x']
    del map_merge_A['uid_y']
    #------ calculate average between expressions lists of different RNAs of the same protein
    group_map_merge_A=map_merge_A.groupby(['ENS_ID'], sort=False).mean(numeric_only=True)
    group_map_merge_A['expression_list_A']= group_map_merge_A.loc[:, group_map_merge_A.columns != 'ENS_ID'].values.tolist()

    map_merge_A=map_merge_A[['uidA', 'uidB','ENS_ID']]
    merge_A=pd.merge(map_merge_A, group_map_merge_A,on='ENS_ID', how='left')
    merge_A=merge_A[['uidA', 'uidB','expression_list_A']]


    #---------------------- B--------------------------#
    map_merge_B= pd.merge(test_mapped, expression_ds_2, left_on='ID_B', right_on='ENS_ID', how='left')
    del map_merge_B['uid_x']
    del map_merge_B['uid_y']
    group_map_merge_B=map_merge_B.groupby(['ENS_ID'], sort=False).mean(numeric_only=True)
    group_map_merge_B['expression_list_B']= group_map_merge_B.loc[:, group_map_merge_B.columns != 'ENS_ID'].values.tolist()

    map_merge_B=map_merge_B[['uidA', 'uidB','ENS_ID']]
    merge_B=pd.merge(map_merge_B, group_map_merge_B,on='ENS_ID', how='left')
    merge_B=merge_B[['uidA', 'uidB','expression_list_B']]


    # DROP DUPLICATES AND! NANs
    #----------------------------------------------------------#
    merge_A = merge_A.drop_duplicates(subset=['uidA','uidB'])
    merge_B = merge_B.drop_duplicates(subset=['uidA','uidB'])
    merge_B =merge_B.drop(['uidA', 'uidB'], axis=1)
    #----------------------------------------------------------#

    concatds=pd.concat([merge_A, merge_B], axis=1)
    exp_list =[]
    for index, rows in concatds.iterrows():
        my_list =[rows.expression_list_A, rows.expression_list_B]
        exp_list.append(my_list)

    calc_list=[]
    for pair in exp_list:
        if (str(pair[0])=='nan')|(str(pair[1])=='nan'):
            calc_list.append('NaN')
        else:
            
            rho, p= spearmanr(pair[0], pair[1])
            calc_list.append(rho)

    return calc_list[0]


irefindex_df = pd.read_csv(DATASETS_DIR / "irefindex_uniprot.csv")
irefindex_df['method'] = irefindex_df['method'].str.extract(r'\((.*?)\)')
def create_uid_dict(df):
    uid_dict = {}
    for idx, row in df.iterrows():
        # Create the key for both (uidA, uidB) and (uidB, uidA) to handle both directions
        key1 = (row['uidA'], row['uidB'])
        key2 = (row['uidB'], row['uidA'])
        value = {'method': row['method'], 'pmids': row['pmids']}
        uid_dict[key1] = value
        uid_dict[key2] = value
    return uid_dict

def check_uids(uid_dict, uidA, uidB):
    # Check for the existence of uidA-uidB or uidB-uidA pairs in the dictionary
    key = (uidA, uidB)
    
    if key in uid_dict:
        exists_value = 1
        method_value = uid_dict[key]['method']
        pmid_value = uid_dict[key]['pmids']
    else:
        exists_value = 0
        method_value = np.nan
        pmid_value = np.nan

    return exists_value, method_value, pmid_value

# Create the dictionary from irefindex
uid_dict = create_uid_dict(irefindex_df)


#-----------------------------------------------------------------------------------------------------------------------#
# Generate and save combinations with features
#-----------------------------------------------------------------------------------------------------------------------#

def generate_combinations_parallel(df, save_dir, max_combinations, num_threads, batch_size, df_to_be_excluded):
    os.makedirs(save_dir, exist_ok=True)

    combinations_seen = set()  # To store seen combinations
    num_combinations = 0

    # Convert df_to_be_excluded to a dictionary for faster lookups
    excluded_pairs = {}
    for _, row in df_to_be_excluded.iterrows():
        uid_a, uid_b = row['uidA'], row['uidB']
        # Store both (uid_a, uid_b) and (uid_b, uid_a) as the pairs are unordered
        excluded_pairs[tuple(sorted([uid_a, uid_b]))] = True

    # Define a function to process a single combination
    def process_combination(combination):
        idx_a, row_a = combination[0]
        idx_b, row_b = combination[1]

        uid_a, uid_b = row_a['uid'], row_b['uid']

        # Sort uids to ensure consistent pair representation
        pair_uids = tuple(sorted([uid_a, uid_b]))

        # Check if pair already seen
        if pair_uids in combinations_seen:
            return None

        # Check if pair exists in FINAL DF TO BE EXCLUDED
        if pair_uids in excluded_pairs:
            return None

        refseq_id_a, refseq_id_b = row_a['refseq_id'], row_b['refseq_id']
        seq_a, seq_b = row_a['seq'], row_b['seq']
        disorder_a = row_a['disorder']
        disorder_b = row_b['disorder']

        # Create disorder_check column
        disorder_check = 1 if disorder_a > 0.5 or disorder_b > 0.5 else 0

        # Get timestamp
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')

        if disorder_check == 1:
            # Generate the row with NaN for feature columns
            nan_values = [np.nan] * 60  # Adjust the number of NaNs based on the number of feature columns
            csv_row = f"{timestamp},{uid_a},{uid_b},{refseq_id_a},{refseq_id_b},{seq_a},{seq_b},{disorder_a},{disorder_b},{disorder_check}," + ','.join(map(str, nan_values)) + "\n"
            csv_row_processed = csv_row  # Processed row is the same in this case
            log_entry = f"{timestamp} - Combination {num_combinations + 1}: {uid_a} - {uid_b} (disorder_check == 1)\n"
        else:
            # Calculate features
            interpro = interpro_calculation(uid_a, uid_b)
            Sequence_similarity = e_value_calc(seq_a, seq_b)
            extra_features = extra_feature_calc(seq_a, seq_b)
            go_sim = go_similarity_code(uid_a, uid_b)
            BP_similarity, MF_similarity, CC_similarity = go_sim
            go_sim_str = ','.join(','.join(map(str, sublist)) for sublist in zip(*go_sim))
            homology_feat = homology_code(refseq_id_a, refseq_id_b, uid_a, uid_b)
            Homologous_in_Mouse, Homologous_in_Drosophila, Homologous_in_Yeast, Homologous_in_Ecoli, Exists_in_MINT, Exists_in_DIP, Exists_in_APID, Exists_in_BIOGRID = homology_feat
            homology_str = ','.join(map(str, homology_feat))
            subcellular = sl_calculation(uid_a, uid_b)
            gd_datasets = gd_datasets_calculation(uid_a, uid_b)
            gd_data_values = gd_datasets[:15]
            gd_data_str = ','.join(map(str, gd_datasets))
            rna1 = rna_dataset_1_calculation(uid_a, uid_b)
            rna2 = rna_dataset_2_calculation(uid_a, uid_b)
            irefindex_result = check_uids(uid_dict, uid_a, uid_b)
            irefindex_check = irefindex_result[0]
            irefindex_method = irefindex_result[1]
            irefindex_pmid = irefindex_result[2]

            # Generate the row to be written to the file
            csv_row = f"{timestamp},{uid_a},{uid_b},{refseq_id_a},{refseq_id_b},{seq_a},{seq_b},{disorder_a},{disorder_b},{disorder_check},{interpro},{Sequence_similarity},{subcellular},{rna1},{rna2},{go_sim_str},{homology_str},{irefindex_check},{irefindex_method},{irefindex_pmid},{gd_data_str},"
            csv_row += ','.join(str(item) for item in extra_features)

            values = [
                timestamp, uid_a, uid_b, refseq_id_a, refseq_id_b, seq_a, seq_b, disorder_a, disorder_b, disorder_check, interpro, Sequence_similarity,
                subcellular, rna1, rna2, BP_similarity, MF_similarity, CC_similarity, Homologous_in_Mouse,
                Homologous_in_Drosophila, Homologous_in_Yeast, Homologous_in_Ecoli, Exists_in_MINT, Exists_in_DIP,
                Exists_in_APID, Exists_in_BIOGRID, irefindex_check, irefindex_method, irefindex_pmid
            ] + gd_data_values + extra_features

            columns = ['timestamp', 'uidA', 'uidB', 'protein_accession_A', 'protein_accession_B', 'seq_A', 'seq_B', 'disorder_a', 'disorder_b', 'disorder_check', 'pfam_interaction', 'Sequence_similarity', 'Subcellular Co-localization?', 'GSE227375_spearman', 'GSE228702_spearman', 'BP_similarity', 'MF_similarity', 'CC_similarity', 'Homologous in Mouse', 'Homologous in Drosophila', 'Homologous in Yeast', 'Homologous in Ecoli', 'Exists in MINT?', 'Exists in DIP?', 'Exists in APID?', 'Exists in BIOGRID?', 'irefindex_check', 'irefindex_method', 'irefindex_pmid', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', 'A %', 'L %', 'F %', 'I %', 'M %', 'V %', 'S %', 'P %', 'T %', 'Y %', 'H %', 'Q %', 'N %', 'K %', 'D %', 'E %', 'C %', 'W %', 'R %', 'G %', 'MW dif', 'Aromaticity dif', 'Instability dif', 'helix_fraction_dif', 'turn_fraction_dif', 'sheet_fraction_dif', 'cys_reduced_dif', 'cys_residues_dif', 'gravy_dif', 'ph7_charge_dif']

            # Create a DataFrame for this combination just for validation
            df_combination = pd.DataFrame([values], columns=columns)

            # Load preprocess models
            loaded_pipeline = joblib.load(DATASETS_DIR / 'models/preprocessing_pipeline.pkl')
            pipeline_features = joblib.load(DATASETS_DIR / 'models/pipeline_features.pkl')

            filtered_new_data = df_combination[pipeline_features]
            labels = df_combination.loc[:, ~df_combination.columns.isin(pipeline_features)]
            for col in filtered_new_data.columns:
                filtered_new_data[col] = filtered_new_data[col].apply(lambda x: x if np.isscalar(x) else np.nan)

            transformed_new_data = loaded_pipeline.transform(filtered_new_data)
            transformed_new_data = pd.DataFrame(transformed_new_data, columns=filtered_new_data.columns)
            final_data = pd.concat([labels, transformed_new_data], axis=1)
            final_data = final_data.fillna(0)

            # Perform predictions
            with open(DATASETS_DIR / "models" / "00001finalSingleModel.pkl.z", "rb") as file:
                model = joblib.load(file)
            filter_file_path = DATASETS_DIR / "models" / "features_selected_1.csv"
            filter_df = pd.read_csv(filter_file_path, index_col=0)
            selected_features = filter_df.columns[filter_df.iloc[0] == 1].tolist()
            df1 = final_data[selected_features]
            predicted_classes = model.predict(df1)
            predicted_probabilities = model.predict_proba(df1)
            class_value = predicted_classes[0]
            classifier_confidence = predicted_probabilities[0][0]

            ####NEW REGRESSION MODELS####ENSEMBLE PREDICTIONS
            def load_model(model_path):
                with open(model_path, 'rb') as file:
                    return joblib.load(file)

            # List of paths to model files and corresponding rows for feature selection
            model_info = [
                {'path': DATASETS_DIR / f'models/Output_dg_dataset_NEW_2/classification_models/000{str(i).zfill(2)}finalSingleModel.pkl.z', 'feature_row': i - 1}
                for i in range(1, 25) if i not in [12] #exclude some models
            ]

            # Load the models
            models = [load_model(info['path']) for info in model_info]

            features_1 = pd.read_csv(DATASETS_DIR / 'models' / 'Output_dg_dataset_NEW_2' / 'feature_selection' / 'features_FinalFront1.csv')
            #print(features_1)

            # Make predictions using the loaded regression models and average them for ensemble
            ensemble_predictions = np.zeros(len(final_data))
            for i, info in enumerate(model_info):
                # Identify columns with '1' in the specified row
                selected_columns = features_1.iloc[info['feature_row']] == 1

                # Select only the columns with '1' in the specified row
                selected_feature_names = features_1.columns[selected_columns]

                # Extract only the selected features from the test dataset
                filtered_test_data = final_data[selected_feature_names]

                # Make predictions using the current model
                predictions = models[i].predict(filtered_test_data)
                ensemble_predictions += predictions

            ensemble_predictions /= len(models)
            regression_value= ensemble_predictions[0]
 
            ########################################

            combinations_seen.add(pair_uids)

            # Update log entry
            log_entry = f"{timestamp} - Combination {num_combinations + 1}: {uid_a} - {uid_b}\n"

            csv_row += f",{class_value},{classifier_confidence},{regression_value}"
            csv_row += '\n'

            # Process values of the first row of filter_data
            filter_values = final_data.iloc[0].values.tolist()
            csv_row_processed = ','.join(str(value) for value in filter_values)
            csv_row_processed += f",{class_value},{classifier_confidence},{regression_value}"
            csv_row_processed += '\n'

        return csv_row, csv_row_processed, log_entry

    def write_results(batch_results):
        with open(os.path.join(save_dir, 'combinations.csv'), 'a') as total_file, \
                open(os.path.join(save_dir, 'combinations_processed.csv'), 'a') as processed_file, \
                open(os.path.join(save_dir, 'combination_log.txt'), 'a') as log_file:
            for csv_row, csv_row_processed, log_entry in batch_results:
                total_file.write(csv_row)
                processed_file.write(csv_row_processed)
                log_file.write(log_entry)

    # Initialize the total and processed files with headers
    with open(os.path.join(save_dir, 'combinations.csv'), 'w') as total_file:
        total_file.write('timestamp,uidA,uidB,protein_accession_A,protein_accession_B,seq_A,seq_B,disorder_a,disorder_b,disorder_check,pfam_interaction,'
                         'Sequence_similarity,Subcellular Co-localization?,GSE227375_spearman,GSE228702_spearman,'
                         'BP_similarity,MF_similarity,CC_similarity,Homologous in Mouse,Homologous in Drosophila,'
                         'Homologous in Yeast,Homologous in Ecoli,Exists in MINT?,Exists in DIP?,Exists in APID?,'
                         'Exists in BIOGRID?,irefindex_check,irefindex_method,irefindex_pmid,0,1,2,3,4,5,6,7,8,9,'
                         '10,11,12,13,14,A %,L %,F %,I %,M %,V %,S %,P %,T %,Y %,H %,Q %,N %,K %,D %,E %,C %,W %,'
                         'R %,G %,MW dif,Aromaticity dif,Instability dif,helix_fraction_dif,turn_fraction_dif, '
                         'sheet_fraction_dif, cys_reduced_dif,cys_residues_dif, gravy_dif,ph7_charge_dif,class_value,classifier_confidence,regression_value\n')

    with open(os.path.join(save_dir, 'combinations_processed.csv'), 'w') as total_file_processed:
        total_file_processed.write('timestamp,uidA,uidB,protein_accession_A,protein_accession_B,seq_A,seq_B,disorder_a,disorder_b,disorder_check,pfam_interaction,'
                                   'Sequence_similarity,Subcellular Co-localization?,GSE227375_spearman,GSE228702_spearman,'
                                   'BP_similarity,MF_similarity,CC_similarity,Homologous in Mouse,Homologous in Drosophila,'
                                   'Homologous in Yeast,Homologous in Ecoli,Exists in MINT?,Exists in DIP?,Exists in APID?,'
                                   'Exists in BIOGRID?,irefindex_check,irefindex_method,irefindex_pmid,0,1,2,3,4,5,6,7,8,9,'
                                   '10,11,12,13,14,A %,L %,F %,I %,M %,V %,S %,P %,T %,Y %,H %,Q %,N %,K %,D %,E %,C %,W %,'
                                   'R %,G %,MW dif,Aromaticity dif,Instability dif,helix_fraction_dif,turn_fraction_dif, '
                                   'sheet_fraction_dif, cys_reduced_dif,cys_residues_dif, gravy_dif,ph7_charge_dif,class_value,classifier_confidence,regression_value\n')

    # Initialize ThreadPoolExecutor with specified number of threads
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        batch = []
        futures = []

        for combination in itertools.combinations(df.iterrows(), 2):
            if num_combinations >= max_combinations:
                break

            futures.append(executor.submit(process_combination, combination))

            if len(futures) >= batch_size:
                for future in as_completed(futures):
                    result = future.result()
                    if result is not None:
                        batch.append(result)

                write_results(batch)
                num_combinations += len(batch)
                batch.clear()
                futures.clear()

        # Final batch
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                batch.append(result)

        write_results(batch)
        num_combinations += len(batch)


def main():
    parser = argparse.ArgumentParser(description="Generate combinations for proteins & Calculate their features")
    parser.add_argument("--input", required=True, help="Path to the input dataset. Must contain 'uid', 'refseq_id', 'seq' & 'disorder' columns. See example input ds")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--max_combinations", type=int, required=True, help="Max number of combinations")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads")
    parser.add_argument("--batch_size", type=int, required=True, help="Batch size for processing")
    parser.add_argument("--additional", help="Path to additional dataset for exclusion (already calculated PPIs)")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"Output directory is set to: {output_path}")

    df = read_and_normalize_input(args.input)

    df_to_be_excluded = create_df_to_be_excluded(additional_df_path=args.additional)
    df_to_be_excluded.to_csv ('df_to_be_excluded.csv', index=False)

    generate_combinations_parallel(
        df=df,
        save_dir=args.output,
        max_combinations=args.max_combinations,
        num_threads=args.threads,
        batch_size=args.batch_size,
        df_to_be_excluded=df_to_be_excluded
    )

if __name__ == "__main__":
    main()