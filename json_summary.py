#!/usr/bin/env python

# This script is intended to load JSON dicts containing resistotypes from Mykrobe predictor,
# It will return a column for each drug, where R = Resistant, S = susceptible, and 0 = unknown.
# The script is an adaptation from this script: https://github.com/iqbal-lab/Mykrobe-predictor/blob/master/scripts/json_to_tsv.py, 
# but returns one line per sample with the resistotypes first and the presence of genes found 
# at least once in the set of samples in subsequent columns.

import argparse
import json
import csv
import os
import sys



parser = argparse.ArgumentParser(description='''Collect and summarize "Mykrobe predict" outputs (JSON files) from a folder into a
single tab separated file (.tsv)''')
parser.add_argument('-i', '--input_dir', type=str,
                    help='Input directory', default="./")
parser.add_argument('-o', '--output_file', type=str,
                    help='Name of the output file', default="summary.tsv")
parser.add_argument('-c', '--min_coverage', type=int,
                    help='Minimum coverage per gene', default=80)
parser.add_argument('-cp', '--copy_number', action="store_true",
                    help='Value written in the output table is "Copy number" instead of coverage')
parser.add_argument('--format', type=str,
                    help='--format', default="wide")

args = parser.parse_args()

json_files = []
for root, dirs, files in os.walk(args.input_dir):
    for name in files:
        if name.endswith(".json"):
            json_files.extend([os.path.join(root, name)])
print (".json-files: ", json_files)




def load_json(f):
    with open(f, 'r') as infile:
        d = json.load(infile)
    return d


def get_drugs(drug_list):
    drugs = []
    for f in args.files:
        try:
            d = load_json(f)
        except ValueError:
            d = {}
        for drug in drug_list:
            if drug not in drugs:
                drugs.append(drug)
    return drugs


def get_phylo_group_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("phylo_group", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_species_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("species", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_lineage_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("lineage", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_file_name(f):
    sample_name = os.path.basename(f).split('.')[0]
    return sample_name


def get_sample_name(f):
    return list(f.keys())[0]

def get_plate_name(f):
    return f.split('/')[-3]

def get_expected_depth(d):
    return str(d.get("expected_depth", -1))


def get_mean_read_length(d):
    return str(d.get("mean_read_length", -1))

def get_called_genes(d, drug=None):
    variants = []
    for variant_name, variant_call in d.items():
        if variant_call.get("_cls") == "Call.SequenceCall":
            per_cov = variant_call.get('info',{}).get('coverage',{}).get("percent_coverage")
            depth = variant_call.get('info',{}).get('coverage',{}).get("median_depth")
            variants.append(":".join([variant_name,
                             str(int(per_cov)),str(int(depth))]))
    return ";".join(variants)

def get_called_genes_wide(d, drug=None):
    variants = {}
    for variant_name, variant_call in d.items():
        variants[variant_name] = {}
        if variant_call.get("_cls") == "Call.SequenceCall":
            per_cov = variant_call.get('info',{}).get('coverage',{}).get("percent_coverage")
            depth = variant_call.get('info',{}).get('coverage',{}).get("median_depth")
            variants[variant_name] ["per_cov"] = per_cov
            variants[variant_name] ["depth"] = depth
    return variants

def get_variant_calls(d):
    variants = []
    for variant_name, variant_call in d.items():
        if variant_call.get("_cls") != "Call.SequenceCall":
            wt_depth = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("median_depth")
            alt_depth = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("median_depth")

            wt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("percent_coverage")
            alt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("percent_coverage")
            if wt_per_cov < 100:
                wt_depth = 0
            if alt_per_cov <100:
                alt_depth =0 
            conf = variant_call.get('info',{}).get('conf',"1")

            variants.append(":".join([variant_name,
                             str(int(alt_depth)),str(int(wt_depth)), str(int(conf))
                           ]))
    return ";".join(variants)

def get_variant_calls_wide(d):
    variants = {}
    for variant_name, variant_call in d.items():
        variants[variant_name] = {}
        if variant_call.get("_cls") != "Call.SequenceCall":
            wt_depth = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("median_depth")
            alt_depth = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("median_depth")
            wt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("reference",{}).get("percent_coverage")
            alt_per_cov = variant_call.get('info',{}).get('coverage',{}).get("alternate",{}).get("percent_coverage")
            if wt_per_cov < 100:
                wt_depth = 0
            if alt_per_cov <100:
                alt_depth =0
            variants[variant_name] ["wt_depth"] = wt_depth
            variants[variant_name] ["alt_depth"] = alt_depth
            variants[variant_name] ["wt_per_cov"] = wt_per_cov
            variants[variant_name] ["alt_per_cov"] = alt_per_cov
            variants[variant_name] ["conf"] = variant_call.get('info',{}).get('conf',"1")
    return variants




if args.format == "long":
    header = [
        "mykrobe_version",
        "file",
        "plate_name",
        "sample",
        "drug",
        "phylo_group",
        "species",
        "lineage",
        "phylo_group_per_covg",
        "species_per_covg",
        "lineage_per_covg", 
        "phylo_group_depth",
        "species_depth",
        "lineage_depth",                  
        "susceptibility",
        "variants (gene:alt_depth:wt_depth:conf)",
        "genes (prot_mut-ref_mut:percent_covg:depth)"]
    with open(os.path.join(args.input_dir, args.output_file), "w+") as output_file:
        output_tsv = csv.writer(output_file, delimiter='\t')
        output_tsv.writerow(header)

        rows = []
        for i, f in enumerate(json_files):
            file = get_file_name(f)
            try:
                d = load_json(f)
                k = list(d.keys())
                d = d[k[0]]
            except ValueError:
                d = {}
            
            phylo_group,phylo_group_per_covg,phylo_group_depth  = get_phylo_group_string(d)
            species,species_per_covg,species_depth  = get_species_string(d)
            lineage,lineage_per_covg,lineage_depth  = get_lineage_string(d)
            sample_name = k[0]
            try:
                plate_name = get_plate_name(f)
            except  IndexError:
                plate_name = ""

            drug_list = sorted(d.get('susceptibility', {}).keys())
            drugs = sorted(drug_list)

            if not drugs:
                drugs = ["NA"]


            for drug in drugs:
                call = d.get('susceptibility', {}).get(drug, {})
                called_by_variants = get_variant_calls(call.get("called_by",{}))
                called_by_genes = get_called_genes(call.get("called_by",{}))
                row = [d.get("version",{}).get("mykrobe-predictor","-1"),
                    file,
                    plate_name,
                    sample_name,
                    drug,
                    phylo_group,
                    species,
                    lineage,
                    phylo_group_per_covg,
                    species_per_covg,
                    lineage_per_covg,                  
                    phylo_group_depth,
                    species_depth,
                    lineage_depth,                
                    call.get(
                        "predict",
                        'N'),
                    called_by_variants, 
                    called_by_genes
                    ]
                # rows.append(row)
                output_tsv.writerow(row)

else:
    drugs_dict = {}
    genes_dict = {}
    variants_dict = {}
        
    for i, f in enumerate(json_files):
        file = get_file_name(f)
        try:
            d = load_json(f)
            k = list(d.keys())
            d = d[k[0]]
        except ValueError:
            d = {}
        drug_list = sorted(d.get('susceptibility', {}).keys())
        drugs = sorted(drug_list)

        if not drugs:
            drugs = ["NA"]
        
        for drug in drugs:
            if not drug in drugs_dict.keys():
                drugs_dict[drug] = ""

            call = d.get('susceptibility', {}).get(drug, {})

            for variant_name, variant_call in call.get("called_by",{}).items():
                if variant_call.get("_cls") != "Call.SequenceCall":
                    if not variant_name in variants_dict.keys():
                        variants_dict[variant_name] = ""

                if variant_call.get("_cls") == "Call.SequenceCall":
                    if not variant_name in genes_dict.keys():
                        genes_dict[variant_name] = ""


    drugs_dict_s = sorted(drugs_dict)
    genes_dict_s = sorted(genes_dict)
    variants_dict_s = sorted(variants_dict)

    column_names = ["file", "sample_name"]

    for key in drugs_dict_s:
        column_names.extend(["S_{}".format(key)])
        # susceptibility_key_order.append(key)
        # susceptibility_key_order = sorted(susceptibility_key_order)

    for key in genes_dict_s:
        column_names.extend(["G_{}".format(key)])
        # seq_calls_key_order.append(key)

    for key in variants_dict_s:
        column_names.extend(["V_{}".format(key)])
        # variant_calls_key_order.append(key)

    with open(os.path.join(args.input_dir, args.output_file), "w+") as output_file:
        output_tsv = csv.writer(output_file, delimiter='\t')
        output_tsv.writerow(column_names)

        for i, f in enumerate(json_files):
            file = get_file_name(f)
            try:
                d = load_json(f)
                k = list(d.keys())
                d = d[k[0]]
            except ValueError:
                d = {}

            sample_name = k[0]
            output_row = [os.path.basename(f), sample_name]

            g = {}
            v = {}
            for key in drugs_dict_s:
                call = d.get('susceptibility', {}).get(key, {})
                output_row.extend(str(call.get("predict", 'N')))
                g[key] = get_called_genes_wide(call.get("called_by",{}))
                v[key] = get_variant_calls_wide(call.get("called_by",{}))

            g_out={}
            for key in genes_dict_s:
                g_out[key] = ""
                for key2, value2 in g.items():
                    # print(value2)
                    if value2.get(key, {}) != {}:
                        a = value2.get(key, {})
                        g_out[key] = str("_".join([str(int(a["per_cov"])),str(int(a["depth"]))]))
            # print(g_out)
            for key, value in sorted(g_out.items()):
                output_row.extend([str(value)])
            
            v_out = {}
            for key in variants_dict_s:
                v_out[key] = ""
                for key2, value2 in v.items():
                    if value2.get(key, {}) != {}:
                        a = value2.get(key, {})
                        v_out[key] = str("_".join([str(int(a["alt_depth"])), str(int(a["wt_depth"])), str(int(a["conf"]))]))

            for key, value in v_out.items():
                output_row.extend([str(value)])
            output_tsv.writerow(output_row)

        legend = "Legend:", "", "", "S_: Susceptibility (S, r, R)", "", "", "", "G_: Tested genes (percent coverage _ depth)", "", "", "", "V_: variants (alt depth _ wt depth _ conf)"
        output_tsv.writerow(legend)
