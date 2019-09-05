from __future__ import print_function

import argparse
import re
import os
from collections import defaultdict, Counter


def _parse_read_name(name):
    # name = name.replace("no_mod", "nomod")
    # name = name.replace("_indel_", "_inside_indel")
    # name = name.replace("_SNP_", "_inside_SNP")
    parts = name.split("_")
    anno = []
    name = name.lower()
    size = 0
    mir = parts[0]
    if name.find("can") > -1:
        anno.append("reference")
        size = 0
    if name.find("ins") > -1:
        anno.append("ins")
        size += 1
    if name.find("del") > -1:
        anno.append("del")
        size += 1
    if name.find("snp") > -1:
        anno.append("snp")
        size += 1
    if name.find("3p-ta") > -1:
        anno.append("3p-ta")
        r = re.search(r"3p-ta[-_][aAtTcGgC]+", name)
        # import pdb; pdb.set_trace()
        size += len(r.group(0).split("-")[-1])
    if name.find("5p-ta") > -1:
        anno.append("5p-ta")
        r = re.search(r"5p-ta[-_][aAtTcGgC]+", name)
        size += len(r.group(0).split("-")[-1])
    if name.find("3p-los") > -1:
        anno.append("3p-los")
        r = re.search(r"3p-los[-_][aAtTgGcC]+", name)
        size += len(r.group(0).split("-")[-1])
    if name.find("5p-los") > -1:
        anno.append("5p-los")
        r = re.search(r"5p-los[-_][aAtTgGcC]+", name)
        size += len(r.group(0).split("-")[-1])
    if name.find("3p-nta") > -1:
        anno.append("3p-nta")
        r = re.search(r"3p-nta[-_][aAtTgGcC]+", name)
        try:
            size += len(r.group(0).split("-")[-1])
        except AttributeError:
            print(name)
            raise ValueError("wrong name: %s" %name)
    if name.find("5p-nta") > -1:
        anno.append("5p-nta")
        r = re.search(r"5p-nta[-_][aAtTgGcC]+", name)
        size += len(r.group(0).split("-")[-1])
    if not anno and not name.startswith("hsa-"):
        anno.append("other")
        size = "NA"
        mir = "other"
    info = {'mir': mir,
            'type': ":".join(anno),
            'size': size}
    if not anno:
        raise ValueError("name is unexpected: %s" % name)
    return info


def _parse_line(line):
    # synthetic = _parse_read_name(line)
    parts = line.split()
    mirna = parts[2]
    # memory[parts[1]] += 1
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
        if line.startswith(">hsa-"):
            # synthetics[line.strip()[1:]] = 0
            seen.add(line.strip()[1:])

print("Total input %s" % len(seen))

summary = defaultdict(list)
with open(args.counts) as inh:
    inh.readline()
    for line in inh:
        data = _parse_line(line)
        _parse_read_name(data[1])
        summary[data[1]].append(data[0])

print("Total counted %s" % len(summary))

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
    results[(auc, detected)]["%s-%s" % (parsed["type"], parsed["size"])] += 1

del summary

for name in seen:
    parsed = _parse_read_name(name)
    results[("FN","None")]["%s-%s" % (parsed["type"], parsed["size"])] += 1

del seen

with open(args.output + "_accuracy.tsv", "w") as outh:
    for (auc, detected) in results:
        for isomir in results[(auc, detected)]:
            print("%s\t%s\t%s\t%s" % (auc, detected, isomir, results[(auc, detected)][isomir]), file=outh)


with open(args.output + "_counts.tsv", "w") as outh:
    for mir in counts_unique:
        print("%s\t%s\t%s" % (mir, counts_unique[mir], counts_all[mir]), file=outh)
