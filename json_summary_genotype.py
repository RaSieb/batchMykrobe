#!/usr/bin/env python
import json
import os
import sys
import csv
import argparse

parser = argparse.ArgumentParser(description='''Collect and summarize "Mykrobe genotype" outputs (JSON files) from a folder into a
single tab separated file (.tsv)''')
parser.add_argument('-i', '--input_dir', type=str,
                    help='Input directory', default="./")
parser.add_argument('-o', '--output_file', type=str,
                    help='Name of the output file', default="summary.tsv")
parser.add_argument('-c', '--min_coverage', type=int,
                    help='Minimum coverage per gene', default=80)
parser.add_argument('-d', '--min_depth', type=int,
                    help='Minimum depth per gene', default=10)
# parser.add_argument('-n', '--copy_number', action="store_true",
#                     help='Value written in the output table is "Copy number" instead of coverage')
parser.add_argument('-v', '--verbous', action="store_true",
                    help='Details on each sample are written in standard output')
# parser.add_argument('-ver', '--version', action="store_true",
#                     help='Value written in the output table is "Version" instead of coverage')

args = parser.parse_args()

def parse_json_in_directory_to_tsv(args):
    json_files = []
    print args
    for root, dirs, files in os.walk(args.input_dir):
        for name in files:
            if name.endswith(".json"):
                json_files.extend([os.path.join(root, name)])
    print json_files

    seq_calls_dict = {}
    for json_file in json_files:
        with open(json_file) as json_data:
            json_dict = json.loads(json_data.read())
            for key, value in json_dict.iteritems():
                sample_cont = value
                if "sequence_calls" in sample_cont:
                    d = json_dict
                    for sample, value in d.iteritems():
                        g = d.get(sample).get("sequence_calls")
                        f = d.get(sample).get("files")
                        if args.verbous:
                            print "+++++", sample, "+++++"
                            print "files:", f
                            print "Sequences with coverage above", args.min_coverage, ":"
                        for seq, seq_info in g.items():
                            seq_info = seq_info[0]
                            if seq_info["info"]["coverage"]["percent_coverage"] >= args.min_coverage & int(seq_info["info"]["coverage"]["median_depth"]) >= args.min_depth:
                                seq_calls_dict[seq] = ""
                                if args.verbous:
                                    if "version" in seq_info["info"]:
                                        print seq, ":", seq_info["info"]["copy_number"], "version:", seq_info["info"]["version"]
                            else:
                                if args.verbous:
                                    print "Cov.", seq_info["info"]["coverage"]["percent_coverage"], "or Depth", seq_info["info"]["coverage"]["median_depth"], "\t(below min!):", seq

    column_names = ["file_name"]
    column_names.extend(["sample_name"])
    seq_calls_key_order = []
    for key in seq_calls_dict:
        column_names.extend(["Sq_{}".format(key)])
        seq_calls_key_order.append(key)
    with open(os.path.join(args.input_dir, args.output_file), "w+") as output_file:
        output_tsv = csv.writer(output_file, delimiter='\t')

        output_tsv.writerow(column_names)
        for json_file in json_files:
            with open(json_file) as json_data:
                output_row = [os.path.basename(json_file)]
                json_dict = json.loads(json_data.read())
                seq_calls_dict = dict.fromkeys(seq_calls_dict, "")

                for sample, value in json_dict.iteritems():
                    output_row.extend([str(sample.split("/")[-1])])
                    if "sequence_calls" in value:
                        for seq, seq_info in json_dict.get(sample).get("sequence_calls").iteritems():
                            seq_info = seq_info[0]
                            if seq_info["info"]["coverage"]["percent_coverage"] >= args.min_coverage:
                                # if args.copy_number:
                                #     seq_calls_dict[seq] = seq_info["info"]["copy_number"]
                                # elif args.version:
                                #     if "version" in seq_info["info"]:
                                #         print seq, ":", seq_info["info"]["copy_number"], "version:", seq_info["info"]["version"]
                                #     seq_calls_dict[seq] = seq_info["info"]["version"]
                                # else:
                                #     seq_calls_dict[seq] = seq_info["info"]["coverage"]["percent_coverage"]
                                cn = seq_info["info"]["copy_number"]
                                if "version" in seq_info["info"]:
                                    v = seq_info["info"]["version"]
                                else: 
                                    v = 0
                                cov = seq_info["info"]["coverage"]["percent_coverage"]
                                d = seq_info["info"]["coverage"]["median_depth"]
                                seq_calls_dict[seq] = str("".join(["c", str(int(cov)), "_d", str(int(d)), "_cn", str(int(cn)), "_v",str(int(v))]))
                    for key in seq_calls_key_order:
                        output_row.extend([str(seq_calls_dict[key])])

                output_tsv.writerow(output_row)
        legend = "Legend:", "", "", "Sq_: gene (variant) tested", "", "", "", "c: percent coverage", "", "", "", "d: depth", "", "", "", "cn: estimated copy number", "", "", "", "v: gene variant"
        output_tsv.writerow(legend)


def main(args):
    parse_json_in_directory_to_tsv(args)

    return


if __name__ == "__main__":
    main(args)
