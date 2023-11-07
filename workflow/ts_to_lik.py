import os, sys
import tskit
import numpy as np
import msprime

from runsmc import liknb
from multiprocessing import Pool
from tqdm import tqdm

TRUE_TREES = '/ceph/users/gbisshop/simplified-trees/n_100/seq_1g'
TRUE_ARGWEAVER_TREES = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-true-trees' 
TSINFER_DIR = '/scratch/gbisshop/output/n_100/seq_1g/tsinfer-trees'
ARGWEAVER_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-simplified-trees'
RELATE_DIR = '/scratch/gbisshop/output/n_100/seq_1g/relate_ts'
TRUE_TREES_OUTPUT = '/ceph/users/gbisshop/output'
TSINFER_OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/tsinfer-trees-likelihoods'
ARGWEAVER_OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-simplified-trees-likelihoods'
RELATE_OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/relate_ts_likelihoods'
ARGWEAVER_TRUE_OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-true-trees-likelihoods'

NE = 10000
MU = 1.25e-8
REC = 1.25e-8

def simplify_ts(ts):
    subsample = list(range(10))
    tss = ts.simplify(samples=subsample, keep_unary=True)
    return simplify_keeping_unary_in_coal(tss)


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


def compute_liks(file, model):
    ts = tskit.load(file)
    if model=='hudson':
        return msprime.log_arg_likelihood(
            ts, 
            REC, 
            NE
        )
    else:
        return liknb.log_likelihood_descending_numba(
            ts, 
            REC, 
            NE
        )


def main():
    app = ['tsinfer', 'relate', 'argweaver']
    app = ['true-argweaver']
    paths = [TSINFER_DIR, RELATE_DIR, ARGWEAVER_DIR]
    paths = [TRUE_TREES]
    out_paths = [TSINFER_OUTPUT_DIR, RELATE_OUTPUT_DIR, ARGWEAVER_OUTPUT_DIR]
    out_paths = [ARGWEAVER_TRUE_OUTPUT_DIR]
    for i, path in enumerate(paths):
        files = os.listdir(path)
        result = np.zeros(len(files))
        idxs = np.zeros(len(files), dtype=int)
        seeds = np.zeros(len(files), dtype=int)
        for j, file in enumerate(files):
            labels = file.split('.')
            seeds[j] = int(labels[0])
            if app[i] == 'argweaver':
                idxs[j] = int(labels[1])
        print(f"[+] Running {app[i]}")
        for j, file in tqdm(enumerate(files), total=len(files)):
            model=None
            if app[i] == 'true-trees':
                model='hudson'      

            result[j] = compute_liks(os.path.join(path, file), model)
        result.dump(os.path.join(out_paths[i], f"{app[i]}.npy"))
        with open(os.path.join(out_paths[i], f'likelihoods-{app[i]}.tsv'), 'w') as file:
            for seed, iteration, value in zip(seeds, idxs, result):
                print(f'{seed}\t{iteration}\t{value}', file=file)


if __name__ == '__main__':
    main()