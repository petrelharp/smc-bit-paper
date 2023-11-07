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

import warnings
warnings.filterwarnings('ignore')

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

TS_DIR = '/scratch/gbisshop/output/n_100/seq_1g/tsinfer-topo'
TRUE_TREES_DIR = '/ceph/users/gbisshop/simplified-trees/n_100/seq_1g' 
TSINFER_OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/tsinfer-topo'
OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/tsinfer-ep-trees-unary'
SAMPLES_DIR = '/scratch/gbisshop/output/n_100/seq_1g/tsinfer-samples'

def run_sim(num_diploids, seq_len, seeds):
    anc_seed, mut_seed = seeds
    ts = msprime.sim_ancestry(
        num_diploids=num_diploids,
        population_size=10**4,
        recombination_rate=1e-8,
        sequence_length=seq_len,
        random_seed=anc_seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=mut_seed)
    ts.dump(TS_DIR + f"-source-{anc_seed}-{mut_seed}.trees")

    # to vcf
    individual_names = [f'HAP{i}' for i in range(ts.num_individuals)]
    with open(f'{anc_seed}-{mut_seed}.vcf', 'w') as output: 
        ts.write_vcf(
            output,
            contig_id='chr1',
            ploidy=2,
            indivual_names=individual_names,
            position_transform=lambda p: [x-1 for x in p],
            )

def run_tsinfer(tsfile):
    ogts = tskit.load(TRUE_TREES_DIR + '/' + tsfile)
    seed = tsfile.rstrip('.trees')
    tsinfer.SampleData.from_tree_sequence(
        ogts, path=SAMPLES_DIR + f"/{seed}.samples", num_flush_threads=2)
    sample_data = tsinfer.load(SAMPLES_DIR + f"/{seed}.samples")
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
    tss.dump(TSINFER_OUTPUT_DIR + f"/{seed}.trees")

def run_tsdate_discrete(tsfile):
    ts = tskit.load(os.path.join(TS_DIR, tsfile))
    seed = tsfile.rstrip('.trees')
    priors = tsdate.build_prior_grid(
        ts, 
        population_size=NE, 
        allow_unary=True, 
        timepoints=1000
    )
    
    dated_ts = tsdate.date(
        ts, 
        priors=priors, 
        mutation_rate=MU
    )   
    # map dates
    dated_ts.dump(OUTPUT_DIR + f"/{seed}.trees")


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


def run_tsdate(tsfile):
    ts = tskit.load(os.path.join(TS_DIR, tsfile))
    seed = tsfile.rstrip('.trees')
    #dated_simpl_tss = tsdate.date(
    #    simpl_tss,
    #    method='variational_gamma',
    #    mutation_rate=MU,
    #    Ne=NE,
    #)
    
    # map dates back
    #new_times = np.zeros(tss.num_nodes)
    #for i, node in enumerate(nodes_map):
    #    if node != -1:
    #        new_times[i] = dated_simpl_tss.nodes_time[node]
    #tables = tss.dump_tables()
    #tables.nodes.time = new_times
    #tables = tss.dump_tables()
    #tables.nodes.time = dated_simpl_tss.nodes_time
    #found = []
    #for edge in tables.edges:
    #    parent_time = dated_simpl_tss.nodes_time[edge.parent]
    #    child_time = dated_simpl_tss.nodes_time[edge.child]
    #    if parent_time <= child_time:
    #        found.append([edge.parent, edge.child, parent_time, child_time])
    #print(found)
    #tables.sort()
    #dated_tss = tables.tree_sequence()
    dated_tss = tsdate.date(
        ts,
        method='variational_gamma',
        mutation_rate=MU,
        Ne=NE,
    )
    dated_tss.dump(OUTPUT_DIR + f"/{seed}.trees")

def main():
    files = os.listdir(TS_DIR)
    with Pool(NUM_CORES) as p:
        for _ in tqdm(p.imap_unordered(run_tsdate, files), total=len(files)):
            pass

if __name__ == "__main__":
    main()