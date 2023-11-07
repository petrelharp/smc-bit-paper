import builtins
import sys
import os
import numpy as np
import json
import msprime
import tskit
import tsinfer
import tsdate

from tqdm import tqdm
from multiprocessing import Pool
from algorithm import *


TS_DIR = ""
SAMPLES_DIR = ""
VCF_DIR = ""

NE=10000 #diploid
MU=1.25e-8
RE=1.25e-8
SEED=1024
ID="chr1"
SEQLEN=1000000
SAMPLES=100
NUM_CORES = 20


TRUE_TS_DIR = '/scratch/gbisshop/output/single/true-ts' 
OUTPUT_DIR = '/scratch/gbisshop/output/single'
SAMPLES_DIR = '/scratch/gbisshop/output/single/single-mut/samples'
TSINFER_OUTPUT_DIR = '/scratch/gbisshop/output/single/single-mut/tsinfer-topo'
TSDATE_OUTPUT_DIR = '/scratch/gbisshop/output/single/single-mut/tsdate-ep'


def run_tsinfer(mts, seed, mut_str):
    tsinfer.SampleData.from_tree_sequence(
        mts, path=OUTPUT_DIR + f"/{mut_str}/samples" + f"/{seed}.samples", num_flush_threads=2)
    sample_data = tsinfer.load(OUTPUT_DIR + f"/{mut_str}/samples" + f"/{seed}.samples")
    ts = tsinfer.infer(sample_data, num_threads=1)
    tables = ts.dump_tables()
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()
    for n in ts.nodes():
        tables.nodes[n.id] = tables.nodes[n.id].replace(
            metadata = json.loads(n.metadata.decode() or "{}"))
    # temporary hack: remove the ultimate ancestor if it exists
    oldest_node = np.argmax(tables.nodes.time)
    if np.sum(tables.edges.parent==oldest_node) == 1:
        # only a single edge connects to the root. This is a unary "ultimate ancestor"
        # and can be removed (it will be removed in later tsinfer versions anyway)
        use = np.arange(tables.nodes.num_rows)
        use = use[use != oldest_node]
        tables.subset(use)
    ts = tables.tree_sequence()
    tss = simplify_keeping_unary_in_coal(ts)
    tss.dump(OUTPUT_DIR + f"/{mut_str}/tsinfer-topo/{seed}.trees")
    return tss

def run_tsdate(ts, seed, mut_str, mut_factor=1):
    tss = ts.simplify()
    dated_tss = tsdate.date(
        tss,
        method='variational_gamma',
        mutation_rate=MU * mut_factor,
        Ne=NE,
    )
    dated_tss.dump(OUTPUT_DIR + f"/{mut_str}/tsdate-ep/{seed}.trees")
    return dated_tss


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

def run_workflow(args):
    ts, label, mut_str, mut_factor = args
    ts_topo = run_tsinfer(ts, label, mut_str)
    ts_dated = run_tsdate(ts_topo, label, mut_str, mut_factor)
    ext_sts, _ = extend_edges(ts_dated)
    ext_sts.dump(OUTPUT_DIR + f"/{mut_str}/tsdate-ep-ext/{label}.trees")    

    return 0

def main():
    files = os.listdir(TRUE_TS_DIR)
    single_ts = files[0]
    ts = tskit.load(os.path.join(TRUE_TS_DIR, single_ts))
    assert single_ts.endswith(".trees")
    label = single_ts.rstrip('.trees')
    rng = np.random.default_rng(9371)
    for mut_str, mut_factor in zip(['single-mut', 'double-mut', 'quadruple-mut'], [1, 2, 4]):
        mut_seeds = rng.integers(1, 2**16, size=100)
        for i in tqdm(range(100)):
            iter_label = label + f'.{i}'
            mts = msprime.sim_mutations(ts, rate=MU*mut_factor, random_seed=mut_seeds[i])
            args = (mts, iter_label, mut_str, mut_factor)
            ret = run_workflow(args)


if __name__ == "__main__":
    main()