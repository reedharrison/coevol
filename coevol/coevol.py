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
	from scipy import linalg
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
		print('Done')

class mca:
	"""Summary

	Attributes
	----------
	X
		data of W binarized by column
	key_col
		List of keys for each column in notation "l.val" where l is the column index from W and val is 
		the unique value that is masked from W[:,l]. Note that l is not unique across keys.
	Z
		data of X transformed
	V
		Matrix of eigenvectors
	N
		Number of rows
	P
		Number of column positions
	J
		Number of non-zero eigenvalues
	f_nS
		Frequency along P
	f_Sp
		Frequency along N
	f_np
		Overall frequency
	y
		Eigenvalues
	y_adj
		Adjusted eigenvalues
	theta
		projection of rows onto standard axes
	psi
		projection of columns onto standard axes
	phi
		mapping of rows to additional feature
	omega
		projection of features onto standard axes


	"""

	def __init__(self, W, phi=None, row=False, col=False):
		"""Summary

		Parameters
		----------
		W
			data set (will be binarized
		L
			true number of 
		phi : (optional)
			feaure matrix describing whether an arbitrary feature is present in each row of X
		proj_row :
			project rows onto principle axes if True
		proj_col :
			project columns onto principle axes if True

		"""
		self.phi = phi
		self.row = row
		self.col = col

		self.L 	= W.shape[1]
		self.N 	= W.shape[0]
		N 		= self.N

		X, self.key_col = binarize(W, lbl=True)
		self.X 			= X
		self.P 			= X.shape[1]
		P 				= self.P

		x_nS = X.sum(axis=0)
		x_Sp = X.sum(axis=1)
		x_SS = X.sum()

		f_nS = x_nS/x_SS
		f_Sp = x_Sp/x_SS
		f_np = X/x_SS

		self.f_nS = f_nS
		self.f_Sp = f_Sp
		self.f_np = f_np

		xv, yv 	= np.meshgrid(f_nS, f_Sp)
		Z 		= f_np / (xv * yv)
		self.Z 	= Z

	def decompose(self):
		"""Summary
		Performs MCA with eigenvalue decomposition

		"""
		Z = self.Z
		L = self.L

		w, v = linalg.eigh(a=np.dot(Z,Z.T))

		chk1 	= w != 1
		chk2 	= w != 0
		chk3 	= w > (1/(L**2))
		mask 	= np.logical_and(chk1, chk2)
		mask 	= np.logical_and(mask, chk3)
		w 		= w[mask]
		v 		= v[:,mask]

		w_adj = ((L/(L-1))**2)*((w-(1/L))**2)

		self.V 		= v
		self.y 		= w
		self.y_adj 	= w_adj
		J 			= len(w)
		self.J 		= J

		# D = np.zeros((J,J))
		# D_adj = np.zeros((J,J))
		# np.fill_diagonal(D, w)
		# np.fill_diagonal(D_adj, w_adj)
		# self.D = D
		# self.D_adj = D_adj

	def percent_variance(self):
		y = self.y_adj
		return 100*y/y.sum()

	def proj_row(self):
		"""Summary
		Projects rows of X onto principle axes

		"""
		N = self.N
		P = self.P
		J = self.J
		V = self.V

		f_nS 	= self.f_nS
		f_Sp 	= self.f_Sp
		y 		= self.y
		y_adj 	= self.y_adj

		R 		= np.zeros((N,N))
		T 		= np.zeros((P,P))
		D 		= np.zeros((J,J))
		D_adj 	= np.zeros((J,J))

		np.fill_diagonal(R, np.sqrt(f_nS))
		np.fill_diagonal(T, np.sqrt(f_Sp))
		np.fill_diagonal(D, np.sqrt(y))
		np.fill_diagonal(D_adj, np.sqrt(y_adj))

		theta 		= np.sqrt(N)*(np.dot(V, D_adj))
		self.theta 	= theta

	def proj_col(self):
		"""Summary
		Projects columns of X onto principle axes

		"""
		N = self.N
		P = self.P
		J = self.J
		Z = self.Z
		V = self.V

		f_nS 	= self.f_nS
		f_Sp 	= self.f_Sp
		y 		= self.y
		y_adj 	= self.y_adj

		R 		= np.zeros((N,N))
		T 		= np.zeros((P,P))
		D 		= np.zeros((J,J))
		D_adj 	= np.zeros((J,J))

		np.fill_diagonal(R, np.sqrt(f_nS))
		np.fill_diagonal(T, np.sqrt(f_Sp))
		np.fill_diagonal(D, np.sqrt(y))
		np.fill_diagonal(D_adj, np.sqrt(y_adj))

		psi 		= np.dot(np.dot(np.dot(T , np.dot(Z.T, V)), D.T), D_adj)
		self.psi 	= psi

	def proj_feat(self):
		"""Summary
		Projects features for rows of X onto principle axes

		"""
		N 	= self.N
		P 	= self.P
		J 	= self.J
		Z 	= self.Z
		V 	= self.V
		phi = self.phi
		if len(phi.shape) == 1:
			C = 1
		if len(phi.shape) > 1:
			C = phi.shape[1]

		f_nS 	= self.f_nS
		f_Sp 	= self.f_Sp
		y 		= self.y
		y_adj 	= self.y_adj

		R 		= np.zeros((N,N))
		H 		= np.zeros((C,C))
		D 		= np.zeros((J,J))
		D_adj 	= np.zeros((J,J))

		np.fill_diagonal(R, np.sqrt(f_nS))
		np.fill_diagonal(C, phi.sum(axis=0))
		np.fill_diagonal(D, np.sqrt(y))
		np.fill_diagonal(D_adj, np.sqrt(y_adj))

		omega = np.dot(np.dot(np.dot(np.dot(np.dot(R, phi), H).T, V), D.T), D_adj)
		self.omega = omega

	def calc(self):
		"""Summary
		Carries out mca analysis. User only need interact with this class method after instantiation of a mca class.

		"""

		self.decompose()
		if self.row:
			self.proj_row()
		if self.col:
			self.proj_col()
		if self.c is not None:
			self.proj_feat()


"""
define coevol functions
"""
def binarize(data, axis=1, lbl=False):
	"""Summary
		Iterate along axis. At each axis position: find all unique values, expand axis element to length of unique values,
			each new element is a mask for a corresponding unique value from the original axis element.

	Parameters
	----------
	data : ndarray
		Numpy array containing data to be binarized
	axis : int(0) or int(1)
		Integer specifying which axis (0=row, 1=column) to binarize along.
	lbl : bool
		When True, this specifies to return the feature labels in a vector

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
	if lbl == True:
		lbls = []
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

		if lbl == True:
			lbls.append(['%s.%s' % (i, x) for x in vals])

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

	if lbl == True:
		lbls = [item for sublist in lbls for item in sublist] # Flatten list
		return (container, lbls)
	else:
		return container