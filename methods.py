#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import pandas as pd
from pre_process import get_pheno, load_npz, merge_geno_pheno, \
						normalization
import os




def main(filenames):
	"""
	"""
	ids, X, snps = [], [], []
	for file in filenames:
		ids_t, X_t, snps_t = load_npz(file)
		ids.append(ids_t)
		X.append(X_t)
		snps.append(snps_t)
	ids = np.concatenate(ids, axis=0)
	X = np.concatenate(X, axis=0)
	merged_df = merge_geno_pheno(ids, X, snps[0])
	print merged_df
	norm_x = normalization(merged_df.ix[:,1:-1])
	print merged_df.columns.values[1:-1]

	

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= 'Several Methods Applied to ..')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--file", dest="file", default=None,
                        help='Pass the path of the file to be process')
    group.add_argument("--folder", dest="folder", default=None,
                        help='Pass the path of the file to be process')
    args = parser.parse_args()
    
    if args.file:
        filenames = [args.file]
    elif args.folder:
        filenames = map(lambda f: os.path.join(args.folder, f), os.listdir(args.folder))
    else:
        parser.print_help()
        sys.exit(1)

    main(filenames)