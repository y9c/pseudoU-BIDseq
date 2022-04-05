import sys
from collections import defaultdict

workdir: config["workdir"]

bed_number = config["bed_number"]
src_dir = config["srcdir"]

data_dir = os.path.dirname(workflow.overwrite_configfiles[-1])
ref_dir = data_dir

if "samples" not in config:
    sys.exit("`samples` is not defined in config file!")

if "references" not in config:
    sys.exit("`references` is not defined in config file!")

REF = config["references"]
# print(workflow.basedir)




group2sample = defaultdict(list)
run2file = {}
sample_ids = []
group_ids = []
run_ids = []
sample2run = defaultdict(list)
read_ids = set()
for g, v in config["samples"].items():
    group_ids.append(g)
    for l in ["input", "treated"]:
        for r, files in v[l].items():
            s = "-".join([g, l, r])
            sample_ids.append(s)
            group2sample[g].append(s)
            for i, r in enumerate(files, 1):
                run_ids.append(s + f"-run{i}")
                sample2run[s].append(s + f"-run{i}")
                run2file[s + f"-run{i}"] = r
                read_ids |= set(r.keys())


rule all:
    input:
        # before combine
