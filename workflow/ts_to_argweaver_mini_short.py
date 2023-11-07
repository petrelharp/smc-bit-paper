# vcf to .sites
import os, sys
import tskit
import numpy as np
import subprocess
import networkx as nx
import pandas as pd
import collections
import msprime
import tsdate
from multiprocessing import Pool

from tqdm import tqdm
from algorithm import *

OUTPUT_DIR = '/scratch/gbisshop/output/argweaver-mini-short'
ARGWEAVER_DIR = '/ceph/users/gbisshop/software/argweaver/bin/arg-sample'
ARGWEAVER_PATH = '/ceph/users/gbisshop/software/argweaver/'
ARGWEAVER_CONV_DIR = '/ceph/users/gbisshop/software/argweaver/bin/smc2arg'

NE = 10000
MU = 2.5e-8
REC = 1.25e-8
NUM_PROCESSES = 1

ANC_SEED = 124511
MUT_SEED = 452121

def run_sim():
    ts = msprime.sim_ancestry(
        samples=10,
        population_size=NE,
        recombination_rate=REC,
        sequence_length=1e4,
        random_seed=ANC_SEED,
        record_full_arg=True
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=MUT_SEED)
    ts.dump(OUTPUT_DIR + f"/true-arg-{ANC_SEED}-{MUT_SEED}.trees")
    return ts


def run_argweaver(sites, max_time):
    ntimes = 100
    seed = sites.rstrip('.sites')
    file = OUTPUT_DIR + f'/{sites}'
    subprocess.run([
        ARGWEAVER_DIR,
        "--sites", file,
        "--output", OUTPUT_DIR + f'/{seed}',
        "--overwrite",
        "--popsize", str(NE),
        "--mutrate", str(MU),
        "--recombrate", str(REC),
        "--iters", "1000",
        "--sample-step", "10",
        "--no-compress-output",
        "--maxtime", f"{max_time}",
        "--ntimes", f"{ntimes}"
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


def simplify_keeping_unary_in_coal(ts, map_nodes=False):
    """
    Keep the unary regions of nodes that are coalescent at least someone in the tree seq
    Temporary hack until https://github.com/tskit-dev/tskit/issues/2127 is addressed
    """
    tables = ts.dump_tables()
    # remove existing individuals. We will reinstate them later
    tables.individuals.clear()
    tables.nodes.individual = np.full_like(tables.nodes.individual, tskit.NULL)

    _, node_map = ts.simplify(map_nodes=True)
    keep_nodes = np.where(node_map != tskit.NULL)[0]
    # Add an individual for each coalescent node, so we can run
    # simplify(keep_unary_in_individuals=True) to leave the unary portions in.
    for u in keep_nodes:
        i = tables.individuals.add_row()
        tables.nodes[u] = tables.nodes[u].replace(individual=i)
    node_map = tables.simplify(keep_unary_in_individuals=True, filter_individuals=False)

    # Reinstate individuals
    tables.individuals.clear()
    for i in ts.individuals():
        tables.individuals.append(i)
    val, inverted_map = np.unique(node_map, return_index=True)
    inverted_map = inverted_map[val != tskit.NULL]
    tables.nodes.individual = ts.tables.nodes.individual[inverted_map]
    if map_nodes:
        return tables.tree_sequence(), node_map
    else:
        return tables.tree_sequence()


def smc_format_to_ts(file):
    name = file.rstrip('.smc')
    subprocess.run([
        "python2",
        ARGWEAVER_CONV_DIR,
        OUTPUT_DIR + '/' + file,
        OUTPUT_DIR + f"/{name}.arg",
    ])
    with open(OUTPUT_DIR + f"/{name}.arg") as f:
        ts = convert_argweaver(f)
        ts = simplify_keeping_unary_in_coal(ts)
        if ts is not None:
            ts.dump(OUTPUT_DIR + f"/{name}-simpl.trees")


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

def modify_argweaver_tree(ts):
    # date with tsdate
    ts = ts.simplify()
    dated_ts = tsdate.date(
        ts,
        method='variational_gamma',
        mutation_rate=MU,
        Ne=NE,
    )
    # extend edges
    ext_sts, _ = extend_edges(dated_ts)
    return ext_sts

def main():
    files = [file for file in os.listdir(OUTPUT_DIR + '/simpl-trees') if file.endswith('simpl.trees')]    
    for file in tqdm(files):
        ts = tskit.load(os.path.join(OUTPUT_DIR + '/simpl-trees', file))
        name = file.rstrip('simpl.trees')
        tss = ts.simplify()
        ext_sts, _ = extend_edges(tss)
        ext_sts.dump(os.path.join(OUTPUT_DIR + '/simpl-ext-trees', f'{name}.simpl.ext.trees'))

    print('---- Finished ----')
    sys.exit()

    for file in tqdm(files):
        # load tree
        ts = tskit.load(os.path.join(OUTPUT_DIR, file))
        name = file.rstrip('simpl.trees')
        new_ts = modify_argweaver_tree(ts)
        new_ts.dump(os.path.join(OUTPUT_DIR, f'{name}-ext.trees'))


if __name__ == '__main__':
    main()
