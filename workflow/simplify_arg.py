import msprime
import sys
import tskit
import os

from tqdm import tqdm

TS_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-args'
OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argweaver-ts-simplified'


def simplify_ts(ts):
	node_type = msprime.NODE_IS_CA_EVENT | msprime.NODE_IS_RE_EVENT
	retain_nodes = np.bitwise_and(ts.nodes_flags, node_type) == 0
	nodes = np.arange(ts.num_nodes)[retain_nodes]
	tss = ts.simplify(
        samples=nodes,
    )
	num_samples = ts.num_samples
	tables = tss.dump_tables()
	# modify tables to correctly mark sample nodes again
	for j in range(num_samples, tss.num_nodes):
	    node_obj = tables.nodes[j]
	    node_obj = node_obj.replace(flags=0)
	    tables.nodes[j] = node_obj
	ts2 = tables.tree_sequence()
	return ts2.simplify(keep_unary=True)


def compute_likelihood(ts):
	pass


def main():
	file_list = [i for i in os.listdir(TS_DIR)]
	for file in tqdm(file_list):
		ts = simplify_ts(tskit.load(TS_DIR + '/' + file))
		ts.dump(OUTPUT_DIR + '/' + file)


if __name__ == "__main__":
	main()