import os, sys
import tskit
import numpy as np
import msprime

from runsmc import liknb
from multiprocessing import Pool
from tqdm import tqdm

SINGLE = '/scratch/gbisshop/output/single/single-mut/'
DOUBLE = '/scratch/gbisshop/output/single/double-mut'
QUADRUPLE = '/scratch/gbisshop/output/single/quadruple-mut'

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


def chunk_ts(ts, windows):
    for window in windows:
        sts = ts.keep_intervals([window], simplify=False).trim()
        sts = simplify_keeping_unary_in_coal(sts)
        yield sts

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

def count_polytomies(ts):
    ret = 0
    for tree in ts.trees():
        ret += np.sum(tree.num_children_array>2)
    return ret

def main():
    paths = [SINGLE, DOUBLE, QUADRUPLE]
    seed = 655987881
    for i, path in enumerate(paths):
        files = os.listdir(path + '/tsdate-ep-ext')
        result = np.zeros(len(files))
        polytomies = np.zeros(len(files), dtype=np.int64)
        idxs = np.zeros(len(files), dtype=int)
        seed = files[0].rstrip().split('.')[0]
        for j, file in tqdm(enumerate(files), total=len(files)):
            labels = file.split('.')
            idxs[j] = int(labels[-2])
            ts = tskit.load(os.path.join(path + '/tsdate-ep-ext', file))
            result[j] = liknb.log_likelihood_descending_numba(
                ts, 
                REC, 
                NE,
                True
            )
            polytomies[j] = count_polytomies(ts)
        result.dump(os.path.join(path + '/liks', f'{seed}_single_liks_rec_corr.npy'))
        with open(os.path.join(path + '/liks', f'{seed}_liks_rec_corr.tsv'), 'w') as file:
            for iteration, value, polytomie in zip(idxs, result, polytomies):
                print(f'{seed}\t{iteration}\t{value}\t{polytomie}', file=file)


if __name__ == '__main__':
    main()