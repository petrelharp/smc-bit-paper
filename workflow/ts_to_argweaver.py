# vcf to .sites
import os, sys
import tskit
import numpy as np
import subprocess
import networkx as nx
import pandas as pd
import collections
import msprime

from multiprocessing import Pool

from tqdm import tqdm

TS_DIR = '/ceph/users/gbisshop/true-trees/n_100/seq_1g'
OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-sites'
ARGWEAVER_SMC_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-smc'
ARGWEAVER_ARG_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-arg'
ARGWEAVER_DIR = '/ceph/users/gbisshop/software/argweaver/bin/arg-sample'
ARGWEAVER_CONV_DIR = '/ceph/users/gbisshop/software/argweaver/bin/smc2arg'
ARGWEAVER_PATH = '/ceph/users/gbisshop/software/argweaver/'
NE = 10000
MU = 1.25e-8
REC = 1.25e-8
NUM_PROCESSES = 20

def run_argweaver(sites):
    seed = sites.rstrip('.sites')
    file = OUTPUT_DIR + f'/{sites}'
    subprocess.run([
        ARGWEAVER_DIR,
        "--sites", file,
        "--output", ARGWEAVER_SMC_DIR + f'/{seed}',
        "--overwrite",
        "--popsize", str(NE),
        "--mutrate", str(MU),
        "--recombrate", str(REC),
        "--iters", "1000",
        "--sample-step", "10",
        "--no-compress-output",
    ])


def convert_argweaver(infile):
    """
    Convert an ARGweaver .arg file to a tree sequence. `infile` should be a filehandle,
    e.g. as returned by the `open` command. An example .arg file is at

    https://github.com/CshlSiepelLab/argweaver/blob/master/test/data/test_trans/0.arg
    """
    start, end = next(infile).strip().split()
    assert start.startswith("start=")
    start = int(start[len("start=") :])
    assert end.startswith("end=")
    end = int(end[len("end=") :])
    # the "name" field can be a string. Force it to be so, in case it is just numbers
    df = pd.read_csv(infile, header=0, sep="\t", dtype={"name": str, "parents": str})
    for col in ("name", "parents", "age"):
        if col not in df.columns:
            raise ValueError(f"Column {col} not found in ARGweaver file")
    name_to_record = {}
    for _, row in df.iterrows():
        row = dict(row)
        name_to_record[row["name"]] = row
    # We could use nx to do this, but we want to be sure the order is correct.
    parent_map = collections.defaultdict(list)

    # Make an nx DiGraph so we can do a topological sort.
    G = nx.DiGraph()
    time_map = {} # argweaver times to allocated time
    for row in name_to_record.values():
        child = row["name"]
        parents = row["parents"]
        time_map[row["age"]] = row["age"]
        G.add_node(child)
        if isinstance(parents, str):
            for parent in row["parents"].split(","):
                G.add_edge(child, parent)
                parent_map[child].append(parent)

    tables = tskit.TableCollection(sequence_length=end)
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()
    breakpoints = np.full(len(G), tables.sequence_length)
    aw_to_tsk_id = {}
    min_time_diff = min(np.diff(sorted(time_map.keys())))
    epsilon = min_time_diff / 1e6
    try:
        for node in nx.lexicographical_topological_sort(G):
            record = name_to_record[node]
            flags = 0
            # Sample nodes are marked as "gene" events
            if record["event"] == "gene":
                flags = tskit.NODE_IS_SAMPLE
                assert record["age"] == 0
                time = record["age"]
            else:
                if record["age"] == 0:
                    time_map[record["age"]] += epsilon
                time = time_map[record["age"]]
                # Argweaver allows age of parent and child to be the same, so we
                # need to add epsilons to enforce parent_age > child_age
                time_map[record["age"]] += epsilon
            tsk_id = tables.nodes.add_row(flags=flags, time=time, metadata=record)
            aw_to_tsk_id[node] = tsk_id
            if record["event"] == "recomb":
                breakpoints[tsk_id] = record["pos"]
    except nx.exception.NetworkXUnfeasible:
        bad_edges = nx.find_cycle(G, orientation="original")
        return None
        #raise nx.exception.NetworkXUnfeasible(
        #    f"Cycle found in ARGweaver graph: {bad_edges}")


    L = tables.sequence_length
    for aw_node in G:
        child = aw_to_tsk_id[aw_node]
        parents = [aw_to_tsk_id[aw_parent] for aw_parent in parent_map[aw_node]]
        if len(parents) == 1:
            tables.edges.add_row(0, L, parents[0], child)
        elif len(parents) == 2:
            # Recombination node.
            # If we wanted a GARG here we'd add an extra node
            x = breakpoints[child]
            tables.edges.add_row(0, x, parents[0], child)
            tables.edges.add_row(x, L, parents[1], child)
            # modify node flags parents
            for parent in parents:
                node_obj = tables.nodes[parent]
                node_obj = node_obj.replace(flags=msprime.NODE_IS_RE_EVENT)
                tables.nodes[parent] = node_obj
        else:
            assert len(parents) == 0
    # print(tables)
    tables.sort()
    # print(tables)
    ts = tables.tree_sequence()
    # The plan here originally was to use the earg_to_garg method to
    # convert the recombination events to two parents (making a
    # standard GARG). However, there are some complexities here so
    # returning the ARG topology as defined for now. There is an
    # argument that we should do this anyway, since that's the structure
    # that was returned and makes very little difference.

    # garg = argutils.earg_to_garg(ts)

    return ts.simplify(keep_unary=True)


def smc_format_to_ts(seed):
    start = 0
    stop = 100
    step = 10
    for i in range(start, stop, step):
        subprocess.run([
            "python2",
            ARGWEAVER_CONV_DIR,
            ARGWEAVER_SMC_DIR + f"/{seed}.{i}.smc",
            ARGWEAVER_ARG_DIR + f"/{seed}.{i}.arg",
        ])
        if os.path.exists(ARGWEAVER_ARG_DIR + f"/{seed}.{i}.arg"):
            with open(ARGWEAVER_ARG_DIR + f"/{seed}.{i}.arg") as f:
                ts = convert_argweaver(f)
                if ts is not None:
                    ts.dump(ARGWEAVER_ARG_DIR + f"/{seed}.{i}.trees")


def subsample_ts(ts, subsample):
    return ts.simplify(samples=subsample)


def to_sites(ts, seed):
    NAMES = [f'n{i}' for i in range(ts.num_samples)]
    with open(OUTPUT_DIR + f"/{seed}.sites", "wt") as file:
        print(
            "NAMES",
            "\t".join(NAMES),
            sep="\t",
            file=file,
        )
        print("REGION",
            "\t".join(["CHR1", str(1), str(int(ts.sequence_length) + 1)]),
            sep="\t",
            file=file,
        )
        for variant in ts.variants():
            print(
                int(variant.site.position),
                "".join(np.array(variant.alleles)[variant.genotypes]),
                sep="\t",
                file=file,
            )


def a():
    files = os.listdir(TS_DIR)
    subsample = list(range(10))
    for file in tqdm(files, total=len(files)):
        ts = tskit.load(TS_DIR + '/' + file)
        ts = subsample_ts(ts, subsample)
        seed = file.rstrip('.trees')
        to_sites(ts, seed)


def b(num_cores):
    files = os.listdir(OUTPUT_DIR)
    print(f"[+] Running argweaver for {len(files)} files.")
    with Pool(num_cores) as p:
        for _ in tqdm(p.imap_unordered(run_argweaver, files), total=len(files)):
            pass


def c(num_cores):
    sys.path.append(ARGWEAVER_PATH)
    files = os.listdir(ARGWEAVER_SMC_DIR)
    seeds = set([i.split('.')[0] for i in files])
    print(f"[+] Running argweaver for {len(seeds)} seeds.")
    with Pool(num_cores) as p:
        for _ in tqdm(p.imap_unordered(smc_format_to_ts, seeds), total=len(seeds)):
            pass


def main():
    c(NUM_PROCESSES)

if __name__ == '__main__':
    main()