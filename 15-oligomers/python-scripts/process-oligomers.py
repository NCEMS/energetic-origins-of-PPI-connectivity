import os, sys
import argparse
import pandas as pd
import numpy as np
import pint
import pint_pandas
from collections import defaultdict

def get_size(row):
    if row['unknown_stoichiometry'] != 1:
        if pd.isna(row['stoichiometry']) or row['stoichiometry'].strip() == '':
            return 0
        values = map(int, row['stoichiometry'].split(';'))
        return sum(values)
    else:
        return np.nan

def ID_homomer(row):
    num_subunits = len(row['complex_string'].split(';'))
    if num_subunits == 1:
        return 1
    else:
        return 0

def ID_homomers_by_size(node, df_complex, N):
    matches = df_complex[df_complex['complex_string'].str.contains(node)]
    for _, row in matches.iterrows():
        if row['oligomer_size'] == N and row['homomer'] == 1:
            return 1
    return 0

def ID_heteromers_by_size(node, df_complex, N):
    matches = df_complex[df_complex['complex_string'].str.contains(node)]
    for _, row in matches.iterrows():
        if row['oligomer_size'] == N and row['homomer'] != 1:
            return 1
    return 0

def build_node_complex_maps(df_complex):
    node_to_complex_ids = defaultdict(set)
    node_to_partners = defaultdict(set)

    for idx, row in df_complex.iterrows():
        if pd.isna(row["complex_string"]):
            continue
        proteins = row["complex_string"].split(";")
        for p in proteins:
            node_to_complex_ids[p].add(idx)
            node_to_partners[p].update(set(proteins) - {p})
    return node_to_complex_ids, node_to_partners

def count_unknown_stoichiometry_complexes(df_complex):
    node_to_unknown_count = defaultdict(int)

    for idx, row in df_complex.iterrows():
        if row.get("unknown_stoichiometry", 0) != 1:
            continue
        if pd.isna(row["complex_string"]):
            continue
        proteins = row["complex_string"].split(";")
        for p in proteins:
            node_to_unknown_count[p] += 1

    return node_to_unknown_count

def is_node_in_any_complex(node, df_complex):
    return int(df_complex['complex_string'].str.contains(node).any())

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--output_dir", default="processed-data")
    parser.add_argument("--organism_tag", default="s288c")
    parser.add_argument("--oligomer_data")
    parser.add_argument("--map_file")
    parser.add_argument("--nodes")
    args = parser.parse_args()

    # read in the nodes_df information
    nodes_df = pd.read_pickle(args.nodes)
    nodes = nodes_df[["node"]].copy()

    # read in the oligomer information from Complex Portal
    olig_df = pd.read_csv(args.oligomer_data, sep="\t")

    # read in and filter the UniProt ID mapping file
    id_map = pd.read_csv(args.map_file, sep="\t", header=None, names=["ID","ID-type","mapped-ID"])
    id_map = id_map[id_map["ID-type"] == "Gene_OrderedLocusName"]

    # process the complex db to get complex participants and stoichiometry for each line in the database
    all_complexes = []
    all_stoich = []
    unknown_stoich = []
    for i, r in olig_df.iterrows():
        protein_list = r["Expanded participant list"].split("|")
        complex_proteins = []
        stoichiometry = []
        unknown = False
        for p in protein_list:
            temp = p.split("(")[1].strip("()")
            complex_proteins.append(p.split("(")[0])
            stoichiometry.append(temp)
            if temp == "0":
                unknown = True
    
        mapped_complex = []
        for p in complex_proteins:
            try:
                mapped = id_map[id_map["ID"] == p]["mapped-ID"].iloc[0]
                mapped_complex.append(mapped)
            except:
                pass

        unknown_stoich.append(unknown)

        #print (i, protein_list, complex_proteins, mapped_complex, stoichiometry, unknown)
        all_complexes.append(mapped_complex)
        all_stoich.append(stoichiometry)

    out = open(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-processed-oligomers-temp.csv", "w")
    out.write("complex_string,stoichiometry,unknown_stoichiometry\n")
    for i in range (0, len(all_complexes)):
        out.write(f"{';'.join(all_complexes[i])},{';'.join(all_stoich[i])},{int(unknown_stoich[i])}\n")
    out.close()

    ###

    # read in the temp data, overwriting the difficult to parse olig_df
    olig_df = pd.read_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-processed-oligomers-temp.csv")

    # compute the size of the oligomer
    olig_df["oligomer_size"] = olig_df.apply(get_size, axis=1)

    # determine whether the oligomer is a homomer or heteromer
    olig_df["homomer"] = olig_df.apply(ID_homomer, axis=1)

    # add information to nodes regarding whether a particular node is involved in different types of interactions

    ## homomers

    # homodimers
    nodes["homodimer"] = nodes["node"].apply(lambda x: ID_homomers_by_size(x, olig_df, 2))

    # homotrimers
    nodes["homotrimer"] = nodes["node"].apply(lambda x: ID_homomers_by_size(x, olig_df, 3))

    # homotetramers
    nodes["homotetramer"] = nodes["node"].apply(lambda x: ID_homomers_by_size(x, olig_df, 4))
    
    ## heteromers

    # heterodimers
    nodes["heterodimer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 2))

    # heterotrimers
    nodes["heterotrimer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 3))

    # heterotetramers
    nodes["heterotetramer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 4))

    # heteropentamers
    nodes["heteropentamer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 5))

    # heterohexamers
    nodes["heterohexamer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 6))

    # heteroheptamers
    nodes["heteroheptamer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 7))

    # heterooctamers
    nodes["heterooctamer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 8))

    # heterononamers
    nodes["heterononamer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 9))

    # heterodecamers
    nodes["heterodecamer"] = nodes["node"].apply(lambda x: ID_heteromers_by_size(x, olig_df, 10))

    ## additional information

    node_to_complex_ids, node_to_partners = build_node_complex_maps(olig_df)
    node_to_unknown_counts = count_unknown_stoichiometry_complexes(olig_df)

    # is the node in any complex
    nodes["in_complex"] = nodes["node"].apply(lambda x: is_node_in_any_complex(x, olig_df))

    # list of nodes this node forms a complex with
    nodes["complex_partners"] = nodes["node"].apply(lambda x: ";".join(sorted(node_to_partners.get(x, set()))))

    # total number of complexes with unknown stoichiometry containing this node
    nodes["complex_count_unknown"] = nodes["node"].apply(lambda x: node_to_unknown_counts.get(x, 0))

    # total number of complexes including this node
    nodes["complex_count"] = nodes["node"].apply(lambda x: len(node_to_complex_ids.get(x, set())))

    ### save the result to file
    nodes.to_csv(f"{args.output_dir}/{args.output_prefix}-{args.organism_tag}-processed-oligomers-per-node.csv", index=False)


if __name__ == "__main__":
    main()
