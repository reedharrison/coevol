print("""
coevol: python module for analysis of protein coevolution - IN DEVELOPMENT!
	version:0.0.0
""")


"""
import python dependencies
"""
try:
	import os as os
	import numpy as np
	import pandas as pd
	from Bio import SeqIO
except:
	print('ERROR: required libraries could not be loaded')


"""
define coevol classes
"""
# class annotateSeqDom:
# 	"""Summary

# 	Attributes
# 	----------
# 	msa_active
# 	msa_family
# 	classes
# 	data

# 	"""

# 	def __init__(self, msa_active, msa_family):#, classes, motifs):
# 		"""Summary

# 		Parameters
# 		----------
# 		msa_active : list
# 			List of Bio.SeqRecord.SeqRecord objects, contains those sequences considered to be active.
# 		msa_family : list
# 			List of Bio.SeqRecord.SeqRecord objects, contains a variety of sequences that may or may not belong to an active class.
# 		classes : list (optional)
# 			List of labels (string or integers) that classify sequences in msa_active (must have same order as appearance of sequences).
# 			If not provided, all active sequences will be labeled uniformly with int(1).
# 		"""

# 	# def cg_annotate(self):
# 	# 	"""Summary
# 	# 	Propogate label from each msa_active sequence to best matching sequence in msa_family.

# 	# 	Parameters
# 	# 	----------
# 	# 	"""
# 	def annotate_by_sequence(self, classes):
# 		"""Summary
# 		Label sequences from msa_active with the corresponding element of classes. Classes must be of the same length of msa_active, or
# 		classes must be of length 1. If classes has only one label, All sequences in msa_active will be labeled uniformly.

# 		Parameters
# 		----------
# 		"""

# 	def set_domains(self, file_domtbl):
# 		"""Summary
# 		Use domtblout from hmmsearch (HH-suite) to identify domains. If annotations have been performed, propogate labels.

# 		Parameters
# 		----------
# 		"""

# 	def annotate_by_domains(self, classes, segments, ref_index=None):
# 		"""Summary
# 		Label all overlapping domains in the sequence alignment with the same classification.


# 		Parameters
# 		----------
# 		"""

# 	def split_domains(self):
# 		"""Summary

# 		Parameters
# 		----------
# 		"""

# 	def build_dom_seqaln(self):
# 		"""Summary

# 		Parameters
# 		----------
# 		"""

# 	def fg_annotate(self):
# 		"""Summary

# 		Parameters
# 		----------
# 		"""

class msa:
	"""Summary

	Attributes
	----------
	seq
	id
	lbl

	Methods
	-------
	mca
	cluster


	"""

	def __init__(self, msa, id=None, lbl=None):
		"""Summary

		Parameters
		----------

		"""

class mca:
	"""Summary

	Attributes
	----------

	"""

	def __init__(self, aln):
		"""Summary

		Parameters
		----------

		"""

"""
define coevol functions
"""
def binarize(data, axis=1):
	"""Summary
		Iterate along axis. At each axis position: find all unique values, expand axis element to length of unique values,
			each new element is a mask for a corresponding unique value from the original axis element.

	Parameters
	----------
	data : ndarray
		Numpy array containing data to be binarized
	axis : int(0) or int(1)
		Integer specifying which axis (0=row, 1=column) to binarize along.

	Tips
	----
	When binarizing a vector, axis parameter won't do anything if the vector shape is (n,). If the vector has been reshaped to something
		similar to (1,n) or (n,1), be sure to specify the axis for binarization.

	"""
	import numpy as np

	shape = data.shape
	if len(shape) == 1:
		n_feat = 1
	else:
		n_feat = data.shape[axis]
	
	container = np.empty(0)
	for i in xrange(n_feat):
		if (n_feat == 1):
			v = data
			axis = np.argmin(v.shape)
		if (axis == 0) and (n_feat != 1):
			v = data[i,:]
		if (axis == 1) and (n_feat != 1):
			v = data[:,i]
		if (axis != 1) and (axis != 0):
			print('ERROR: the axis parameter must be 0 or 1')

		vals = np.unique(v)
		n_val = len(vals)
		n_rec = len(v)

		if (axis == 0):
			subset = np.zeros((n_val, n_rec))
			for j in xrange(n_val):
				subset[j,np.where(v == vals[j])] = 1

		if (axis == 1):
			subset = np.zeros((n_rec, n_val))
			for j in xrange(n_val):
				subset[np.where(v == vals[j]), j] = 1

		if container.size != 0:
			container = np.concatenate((container, subset), axis=axis)
		if container.size == 0:
			container = subset

	return container