#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
from sklearn.svm import SVC
from sklearn.decomposition import RandomizedPCA
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.cross_validation import train_test_split
from pre_process import load_npz, merge_geno_pheno, RawDataPheno, \
						normalization
from plot_tools import plot_PCA, plot_confusion_matrix
from sklearn.metrics import accuracy_score, classification_report
from pystruct.models import GraphCRF
from pystruct.learners import OneSlackSSVM


def get_PCA(X):
	"""
	"""
	# Project the data to a 2D space for visualization
	Xp = RandomizedPCA(n_components=2, random_state=1).fit_transform(X)
	return Xp

def conditional_random_fields(X, y):
	"""
	"""

	X_ = [(np.atleast_2d(x), np.empty((0, 2), dtype=np.int)) for x in X]
	Y = y.reshape(-1, 1)

	X_train, X_test, y_train, y_test = train_test_split(X_, Y)

	pbl = GraphCRF()
	svm = OneSlackSSVM(pbl)

	svm.fit(X_train, y_train)
	y_pred = np.vstack(svm.predict(X_test))
	print("Score with pystruct crf svm: %f "
	      % (np.mean(y_pred == y_test)))
	print classification_report(y_test, y_pred)
	plot_confusion_matrix(y_test, y_pred)

def extra_tree(X_train, y_train):
	"""
	"""
	clf = ExtraTreesClassifier(n_estimators=4,
							   criterion='gini',
	                           max_features=0.2, 
	                           n_jobs=4,
	                           random_state=1).fit(X_train, y_train)
	#y_pred = clf.predict(X_test)

	return clf

def support_vector(X_train, y_train):
	"""
	"""
	clf = SVC(kernel='linear')
	clf.fit(X_train, y_train)
	
	return clf

def main(filenames, filename_pheno, phenos):
	"""
	"""
	ids, X = [], []
	snps = None
	for file in filenames:
		ids_t, X_t, snps = load_npz(file)
		ids.append(ids_t)
		X.append(X_t)
	# This should go away


	ids = np.concatenate(ids, axis=0)
	X = np.concatenate(X, axis=0)

	data_pheno_ear = RawDataPheno(filename_pheno, phenos)
	df_pheno_ear = data_pheno_ear.get_pheno()

	
	merged_df = merge_geno_pheno(df_pheno_ear, ids, X, snps)
	#print merged_df
	norm_x = normalization(merged_df.ix[:,2:])
	print merged_df.ix[:,2:]
	
	X_pca = get_PCA(norm_x)
	#merged_df[data_pheno.pheno_name[1]].replace(to_replace=-1, value=0, inplace=True)
	merged_df[data_pheno_ear.pheno_name[-1]] = merged_df[data_pheno_ear.pheno_name[-1]].astype(np.int)
	y = merged_df[data_pheno_ear.pheno_name[-1]].values
	print "Pheno working : {} on {}".format(data_pheno_ear.pheno_name[-1], merged_df.columns)
	"""if np.nan in y or -1 in y:
		print "hay nan :("
	else:
		print "la muestra esta completa"
	"""
	plot_PCA(X_pca, y)
	
	#conditional_random_fields(norm_x, y)
	X_train, X_test, y_train, y_test = train_test_split(norm_x, y)
	
	clf = support_vector(X_train, y_train)
	clf = extra_tree(X_train, y_train)
	
	y_pred = clf.predict(X_test)
	clf.score(X_test, y_test)
	print classification_report(y_test, y_pred)
	print data_pheno_ear.pheno_name[-1]
	plot_confusion_matrix(y_test, y_pred)


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
    
    #TODO pass by parameter

    filename_pheno = 'data/prediction/edad_sexo_oreja.txt'
    #phenos = ['EarProtrusion', 'LobeAttachment', 'LobeSize', 'HelixRolling',
    #		  'Crushelixexpression', 'SuperiorCrusofantihelixexpression',
    #		  'Foldofantihelix', 'Darwinstubercle', 'Tragussize', 'Antitragussize']
    phenos = ['LobeAttachment', 'Crushelixexpression', 'Tragussize']
    for pheno in phenos:
    	phenos = ['IID', 'sex', pheno]
    	main(filenames, filename_pheno, phenos)