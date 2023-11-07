import os, sys
import tskit
import numpy as np
import msprime

from runsmc import liknb
from multiprocessing import Pool
from tqdm import tqdm

BASE_DIR_MINI = '/scratch/gbisshop/output/argweaver-mini'
BASE_DIR_MINI_SHORT = '/scratch/gbisshop/output/argweaver-mini-short'


NE = 10000
MU = 1.25e-8
REC = 1.25e-8


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


def compute_liks(ts, model):
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
    filedirs = [BASE_DIR_MINI, BASE_DIR_MINI_SHORT]
    subdirs = ['/simpl-trees', '/simpl-ext-trees', '/ext-trees']

    for filedir in filedirs:
        files = os.listdir(filedir)
        true_ts = None
        for file in files:
            if file.startswith('true-arg'):
                true_ts = file
                break
        ts = tskit.load(os.path.join(filedir, true_ts))
        tss = simplify_keeping_unary_in_coal(ts)
        liks = [compute_liks(tss, 'smc'), compute_liks(ts, 'hudson')]
        with open(os.path.join(filedir, f'true-liks.tsv'), 'w') as file:
            print(f'smc\t{liks[0]}', file=file)
            print(f'hudson\t{liks[1]}', file=file)

        for subdir in subdirs:    
            files = [file for file in os.listdir(filedir + subdir) if file.endswith('.trees')]
            result = np.zeros(len(files))
            idxs = np.zeros(len(files), dtype=int)
            seeds = np.zeros(len(files), dtype=int)
            for j, file in enumerate(files):
                labels = file.split('.')
                seeds[j] = int(labels[0])
                temp_idx = ''.join([x for x in labels[1] if x.isnumeric()])
                idxs[j] = int(temp_idx)
            for j, file in tqdm(enumerate(files), total=len(files)):
                ts = tskit.load(os.path.join(filedir + subdir, file))
                result[j] = compute_liks(ts, 'smc')
            result.dump(os.path.join(filedir + subdir, f"liks.npy"))
            with open(os.path.join(filedir + subdir, f'liks.tsv'), 'w') as file:
                for seed, iteration, value in zip(seeds, idxs, result):
                    print(f'{seed}\t{iteration}\t{value}', file=file)


if __name__ == '__main__':
    main()