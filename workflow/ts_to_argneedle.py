# vcf to .sites
import os
import tskit
import numpy as np
import subprocess

from tqdm import tqdm

INPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/temp'
INFO_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argneedle-info'
OUTPUT_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argneedle-arg'
ARGNEEDLE_TS_DIR = '/scratch/gbisshop/output/n_100/seq_1g/argneedle-trees'
TS_DIR = '/ceph/users/gbisshop/true-trees/n_100/seq_1g'

NE = 10000
MU = 1.25e-8
REC = 1.25e-8
NUM_PROCESSES = 20

def run_argneedle(seed):
    subprocess.run([
        "arg_needle",
        "--hap_gz", INPUT_DIR + f'/{seed}.haps',
        "--map", INPUT_DIR + f'/{seed}.map',
        "--normalize", "1",
        "--normalize_demography", INFO_DIR + '/demo.demo',
        "--out", OUTPUT_DIR + f'/{seed}'
    ])


def generate_demo():
    with open(INFO_DIR + '/demo.demo', 'w') as file:
        print(f'0.0\t{NE * 2}')

def generate_map(seed):
    ts = tskit.load(TS_DIR + f'/{seed}.trees')
    with open(INFO_DIR + f'/{seed}.map', 'w') as file:
        for variant in ts.variants():
            if variant.num_alleles <= 2:
                # chr1, snp_id, cM, bp
                rec_pos = variant.position * (REC * 100)
                print(f'1\tsnp{variant.site.id}\t{rec_pos}\t{variant.position+1}', file=file)

def convert_argweaver(infile):
    pass


def main():
    seeds = [file.rstrip('.trees') for file in os.listdir(TS_DIR)]
    for seed in tqdm(seeds, total=len(seeds)):
        run_argneedle(seed)


if __name__ == '__main__':
    main()