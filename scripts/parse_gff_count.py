from __future__ import print_function

import argparse
import os
from collections import defaultdict, Counter


def _parse_read_name(name):
    name = name.replace("no_mod", "nomod")
    name = name.replace("_indel_", "_inside_indel")
    name = name.replace("_SNP_", "_inside_SNP")
    parts = name.split("_")
    anno = []
    name = name.lower()
    if name.find("nomod") > -1:
        anno.append("reference")
    if name.find("indel") > -1:
        anno.append("indel")
    if name.find("SNP") > -1:
        anno.append("SNP")
    if name.find("3prime_templated") > -1:
        anno.append("3prime_templated")
    if name.find("5prime_templated") > -1:
        anno.append("5prime_templated")
    if name.find("3prime_loss") > -1:
        anno.append("3prime_loss")
    if name.find("5prime_loss") > -1:
        anno.append("5prime_loss")
    if name.find("3prime_addition") > -1:
        anno.append("3prime_non")
    if name.find("5prime_addition") > -1:
        anno.append("5prime_non")
    if name.find("3prime_non") > -1:
        anno.append("3prime_non")
    if name.find("5prime_non") > -1:
        anno.append("5prime_non")
    info = {'mir': parts[0],
            'type': ":".join(anno)}
    return info


def _parse_line(line, memory):
    # synthetic = _parse_read_name(line)
    parts = line.split()
    mirna = parts[2]
    memory[parts[1]] += 1
    # return {'mirna': mirna, 'name': parts[1]}
    return [mirna, parts[1]]


parser = argparse.ArgumentParser()
parser.add_argument("--counts",
                    help="File with gff count data.", required=True)
parser.add_argument("--input",
                    help="File with input fasta file.", required=True)
parser.add_argument("--output",
                    help="prefix for output files.", required=True)

args = parser.parse_args()
synthetics = Counter()
counts_unique = Counter()
counts_all = Counter()
seen = set()

with open(args.input) as inh:
    for line in inh:
        if line.startswith(">"):
            synthetics[line.strip()[1:]] = 0
            seen.add(line.strip()[1:])

summary = defaultdict(list)
with open(args.counts) as inh:
    for line in inh:
        data = _parse_line(line, synthetics)
        _parse_read_name(data[1])
        summary[data[1]].append(data[0])

results = defaultdict(Counter)
for read in  summary:
    if read in seen:
        seen.remove(read)
    for hit in summary[read]:
        parsed = _parse_read_name(read)
        if parsed["mir"] == hit:
            auc = "TP"
            break
        else:
            #if parsed["type"] == "reference":
            #    print("reference non detected %s-%s-" % (name, parsed["type"], hit, parsed["mir"]))
            auc = "FP"
    if len(summary[read]) == 1:
        detected = "once"
        counts_unique[hit] += 1
        counts_all[hit] += 1
    else:
        counts_all[hit] += 1
        detected = "multiple"
    results[(auc, detected)][parsed["type"]] += 1

del summary

for name in seen:
    parsed = _parse_read_name(name)
    results[("FN","None")][parsed["type"]] += 1

del seen

with open(args.output + "_accuracy.tsv", "w") as outh:
    for (auc, detected) in results:
        for isomir in results[(auc, detected)]:
            print("%s\t%s\t%s\t%s" % (auc, detected, isomir, results[(auc, detected)][isomir]), file=outh)


with open(args.output + "_counts.tsv", "w") as outh:
    for mir in counts_unique:
        print("%s\t%s\t%s" % (mir, counts_unique[mir], counts_all[mir]), file=outh)
