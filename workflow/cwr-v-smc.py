import msprime
import sys
import tskit

NE="10000" #diploid
MU="1.25e-8"
RE="1.25e-8"
SEED="1024"
ID="chr1"
SEQLEN=1000000
SAMPLES=100

OUTPUT_DIR = '../true-trees'

def sim_ts(seed):
	ts = msprime.sim_ancestry(
		samples=SAMPLES,
		population_size=NE,
		random_seed=seed,
		model='hudson',
		record_full_arg=True,
		sequence_length=SEQLEN
	)
	return ts

def simplify_ts(ts):
	node_type = msprime.NODE_IS_CA_EVENT | msprime.NODE_IS_RE_EVENT
	retain_nodes = np.bitwise_and(ts.tables.nodes.flags, node_type) == 0
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
	if len(sys.argv) == 4:
		unique_id = sys.argv[1]
		num_replicates = sys.argv[2]
		start_seed = sys.argv[3]
	else:
		unique_id = "test_sim"
		start_seed = 42
		num_replicates = 10
	
	# Set a fixed random seed for reproducibility
	rng = np.random.default_rng(start_seed)
	random_seeds = rng.integers(1, 2**31, size=num_replicates)
    results = np.zeros((2, random_seeds.size), dtype=np.float64)

    for i in random_seeds.size:
    	ts = sim_ts(random_seeds[i])
    	results[0, i] = msprime.log_arg_likelihood(
    		ts, 
    		Ne=NE,
    		recombination_rate=RE,
    	)
    	tss = simplify_ts(ts)
    	results[1,i] = compute_likelihood(tss)

    output = f'results_{unique_id}_{num_replcates}_{start_seed}.tsv'
    with open(output, 'a') as outfile:
    	for seed, value1, value2 in zip(random_seeds, results[0], results[1]):
    		print(f'{seed}\t{value1}\t{value2}', file=outfile)

if __name__ == "__main__":
	main()