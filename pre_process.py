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


class RawDataPheno(object):
    """
    docstring for RawDataPheno
    """
    def __init__(self, filename, pheno_name):
        super(RawDataPheno, self).__init__()
        self.filename = filename
        self.pheno_name = pheno_name
    
    #@profile
    def get_pheno(self):
        """
        docstring for get_pheno
        """
        df = pd.read_csv(self.filename, usecols=self.pheno_name, sep=r"\t")
        return df   

class RawDataGeno(object):
    """
    docstring for RawData
    """
    def __init__(self, filenames, snps_to_use):
        """
        """
        super(RawDataGeno, self).__init__()
        self.filenames = filenames
        # Get the header from the first file
        self.all_snps = self.get_all_snps(filenames[0])
        self.snps_to_use = self.get_filter_snps(snps_to_use)
        self.snps_to_use = [snp[0] for snp in self.snps_to_use]
    

    #@profile
    def get_all_snps(self, filename):
        """
        Return all the columns of the original file.
        """

        df = pd.read_csv(filename, sep=r"\s*", nrows=1)
        cols = df.columns.values
        del df
        return cols
    
    #@profile
    def get_filter_snps(self, snps_to_use):
        """
        Get only the columns for the ID and wanted snps.
        """
        #TODO get by parameters the snps to filter
        union_term = "_\w|"
        search_term = "(IID|" + union_term.join(snps_to_use) + "_\w)"
        search_re = re.compile(search_term).search
    
        return [ ( l, m.group(1) ) for l in self.all_snps for m in (search_re(l),) if m]
    
    #@profile
    def get_geno(self, filename):
        """
        Return only the values of the selected snps and the
        ID of the individual.
        """
        df = pd.read_csv(filename, header=None, names=self.all_snps,
                          usecols=self.snps_to_use, sep=r"\s*")
        
        return df

    #@profile
    def save_binary_file(self, df, file):
        """
        Take geno values, snps and id to npz files
        """
        #TODO check if this have to be a method

        X = df.ix[:, 1:].values
        snps = self.snps_to_use
        ids = df['IID'].values
        np.savez("out_".format(os.path.splitext(os.path.basename(file))[0]), 
                               X=X, snps=snps, ids=ids)


###### CLEAN THIS #######
#@profile
def merge_geno_pheno(df_pheno, ids, X, snps):
    """
    """
    # TODO pass by parameters

    #df_pheno = get_pheno('data/prediction/EarPhenoNew2NA.txt', ['IID', 'LobeSize'])
    df_geno_id = pd.DataFrame(ids, columns=['IID'])
    df_geno_data = pd.DataFrame(X, columns=snps[1:])
    
    df_geno = pd.concat([df_geno_id, df_geno_data], axis=1)
    df_merged = pd.merge(df_geno, df_pheno, on='IID')
    df_merged.fillna(-1, inplace=True)
    return df_merged

#@profile
def load_npz(filename):
    """
    """
    data = np.load(filename)
    ids = data['ids']
    X = data['X']
    snps = data['snps']

    return ids, X, snps

def normalization(X):
    """
    """    
    imputer = Imputer(missing_values=-1, strategy="most_frequent")
    X = imputer.fit_transform(X)

    return X


###### MAIN PART #######
#@profile
def main(filenames, snps_to_use):
    """
    Save the raw data into npz files for faster manipulation
    """
    #TODO file by parameter
    my_data = RawDataGeno(filenames, snps_to_use)
    for file in filenames:
        df = my_data.get_geno(file)
        my_data.save_binary_file(df, file)

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
        filenames = sorted(map(lambda f: os.path.join(args.folder, f), os.listdir(args.folder)))
    else:
        parser.print_help()
        sys.exit(1)
    #TOD pass this by parameter
    snps_to_use = ["rs17023457", "rs3827760", "rs2080401", "rs10212419", 
                   "rs1960918", "rs263156", "rs1619249"]
    main(filenames, snps_to_use)