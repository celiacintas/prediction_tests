#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
import argparse
import os
import re

# If you have bed bim fam files first run
# plink --bfile bed_file --recodeA --out raw_file
# this return the ped and map files

#@profile
def merge_geno_pheno(ids, X, snps):
    """
    """
    # TODO pass by parameters

    df_pheno = get_pheno('data/prediction/EarPhenoNew2NA.txt', ['IID', 'LobeSize'])
    df_geno_id = pd.DataFrame(ids, columns=['IID'])
    df_geno_data = pd.DataFrame(X, columns=snps[1:])
    
    df_geno = pd.concat([df_geno_id, df_geno_data], axis=1)
    df_merged = pd.merge(df_geno, df_pheno, on='IID')
    df_merged.fillna(-1, inplace=True)

    return df_merged

#@profile
def get_pheno(filename, filter_columns):
    """
    """
    df = pd.read_csv(filename, usecols=filter_columns, sep=r"\t")
    return df

#@profile
def load_npz(filename):
    """
    """
    data = np.load(filename)
    ids = data['ids']
    X = data['X']
    snps = data['snps']

    return ids, X, snps

#@profile
def get_geno(filename, header, filter_columns):
    """
    Return only the values of the selected snps and the
    ID of the individual.
    """
    df = pd.read_csv(filename, header=None, names=header, usecols=filter_columns, sep=r"\s*")
    return df

#@profile
def get_filter_snps(list_all_snps):
    """
    Get only the columns for the ID and wanted snps.
    """
    #TODO get by parameters the snps to filter
    search_re = re.compile('(IID|rs17023457_\w|rs3827760_\w|rs2080401_\w|rs10212419_\w|rs1960918_\w|rs263156_\w |rs1619249_\w)').search
    
    return [ ( l, m.group(1) ) for l in list_all_snps for m in (search_re(l),) if m]


#@profile
def get_all_snps(filename):
    """
    Return all the columns of the original file.
    """

    df = pd.read_csv(filename, sep=r"\s*", nrows=1)
    cols = df.columns.values
    del df
    return cols

#@profile
def save_binary_file(df, file):
    """
    Take geno values, snps and id to npz files
    """

    X = df.ix[:, 1:].values
    snps = df.columns.values
    ids = df['IID'].values
    np.savez("data/kaustubh-all-numpy_{}.npz".format(os.path.splitext(os.path.basename(file))[0]), X=X, snps=snps, ids=ids)
    print "save numpy binary"


def normalization(X):
    """
    """    
    print "Before imputation:", np.unique(X)
 
    imputer = Imputer(missing_values=-1, strategy="most_frequent")
    X = imputer.fit_transform(X)

    print "After imputation:", np.unique(X), X.shape

    return X

#@profile
def main(filenames):
    """
    Save the raw data into npz files for faster manipulation
    """
    #TODO file by parameter
    header = get_all_snps('data/prediction/split_rawaa')
    filter_snps = get_filter_snps(header)
    filter_snps = [snp[0] for snp in filter_snps]

    for file in filenames:
        df = get_geno(file, header, filter_snps)
        save_binary_file(df, file)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
                                     'Pre Processing Data')
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