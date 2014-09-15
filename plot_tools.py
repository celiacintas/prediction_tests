#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report

def plot_PCA(Xp, y):
	"""
	"""
	# Plot individuals 
	populations = np.unique(y)

	colors = plt.get_cmap("hsv")
	plt.figure(figsize=(10, 4))
	hair = np.unique(y)
	for i, p in enumerate(hair):
	    mask = (y == p)
	    plt.scatter(Xp[mask, 0], Xp[mask, 1], 
	                c=colors(1. * i / 11), label=p)
	    
	#plt.xlim([-30, 50])
	plt.legend(loc="best")
	plt.show()


def plot_confusion_matrix(y_pred, y):
	"""
	"""
	plt.imshow(confusion_matrix(y, y_pred), cmap=plt.cm.binary, interpolation='nearest')
	plt.colorbar()
	plt.xlabel('true value')
	plt.ylabel('predicted value')
	plt.show()