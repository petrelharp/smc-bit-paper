import os, sys
import tskit
import numpy as np
import msprime
import pyslim

from tqdm import tqdm

SLIM_TS = '/scratch/gbisshop/bs'

NE = 10000
MU = 2.5e-6
REC = 1.25e-8
S = -0.0025
L = 1e6
DELTA_POS = 1E3


def simplify_keeping_unary_in_coal(ts, map_nodes=False):
    """
    Keep the unary regions of nodes that are coalescent at least someone in the tree seq
    Temporary hack until https://github.com/tskit-dev/tskit/issues/2127 is addressed
    """
    tables = ts.dump_tables()
    # remove existing individuals. We will reinstate them later
    tables.individuals.clear()
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
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


def main():
    rng = np.random.default_rng(101)
    random_seed = rng.integers(0, 2**16)
    ts = tskit.load(os.path.join(SLIM_TS, 'bs.co.false.trees'))
    random_samples = rng.choice(np.arange(20000), replace=False, size=200)
    
    tss = ts.simplify(random_samples, keep_unary=True)
    cotss = simplify_keeping_unary_in_coal(tss)
    cotss.dump(os.path.join(SLIM_TS, 'bs.co.false.subsample.trees'))

    height_for_pos = np.zeros(int(cotss.sequence_length))
    for tree in cotss.trees():
        mean_height = np.mean([tree.time(root) for root in tree.roots])
        left, right = map(int, tree.interval)
        height_for_pos[left: right] = mean_height
    height_for_pos.dump('tree_height_bs.npy')


if __name__ == '__main__':
    main()