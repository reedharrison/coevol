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
	from scipy import linalg, stats
	from scipy.cluster.vq import kmeans2
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
	key_lbl
		List of keys for sequence labels in same format as key_col
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

	def __init__(self, W, phi=None, row=False, col=False, tol=0.00001):
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
		tol:
			eigenvalues falling within the range of 1-tol<y<1+tol will be discarded as trivial

		"""
		self.phi, self.key_lbl = binarize(phi, axis=1, lbl=True)
		self.row = row
		self.col = col
		self.tol = tol

		self.L 	= W.shape[1]
		self.N 	= W.shape[0]
		N 		= self.N

		X, self.key_col = binarize(W, lbl=True)
		self.X 			= X
		self.P 			= X.shape[1]
		P 				= self.P

		x_nS = X.sum(axis=1) 	# Row sums
		x_Sp = X.sum(axis=0) 	# Col sums
		x_SS = X.sum()			# Array sum

		f_nS = x_nS/x_SS
		f_Sp = x_Sp/x_SS
		f_np = X/x_SS

		self.f_nS = f_nS
		self.f_Sp = f_Sp
		self.f_np = f_np

		xv, yv 	= np.meshgrid(f_Sp, f_nS)
		Z 		= f_np / np.sqrt(xv * yv)
		self.Z 	= Z

	def decompose(self):
		"""Summary
		Performs MCA with eigenvalue decomposition

		"""
		Z = self.Z
		L = self.L
		tol = self.tol

		w, v = linalg.eigh(a=np.dot(Z,Z.T))

# Notes: filtering out trivial eigenvalues was not working without tolerance implementation
		# chk1 	= w != 1
		# chk2 	= w != 0
		# chk3 	= w > (1/(L**2))
		# mask 	= np.logical_and(chk1, chk2)
		# mask 	= np.logical_and(mask, chk3)
		# w 		= w[mask]
		# v 		= v[:,mask]

		chk1 	= np.logical_or(w < 1-tol, w > 1+tol)
		chk2 	= w >= (1/(L**2))
		mask 	= np.logical_and(chk1, chk2)
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

	def percent_variance(self, mode=None):
		y = self.y_adj
		variance = 100*y/y.sum()
		if mode is None:
			return np.flipud(variance)
		else:
			return np.flipud(variance)[mode]

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

		R 		= np.mat(R)
		T 		= np.mat(T)
		D 		= np.mat(D)
		D_adj 	= np.mat(D_adj)
		V 		= np.mat(V)

		# theta 		= np.sqrt(N)*(np.dot(V, D_adj))
		theta 		= np.sqrt(N)*(V * D_adj)
		theta_std	= np.sqrt(N)*(V * D.I * D_adj)
		self.theta 	= theta.A
		self.theta_std = theta_std.A

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

		R 		= np.mat(R)
		T 		= np.mat(T)
		D 		= np.mat(D)
		D_adj 	= np.mat(D_adj)
		Z 		= np.mat(Z)
		V 		= np.mat(V)

		np.fill_diagonal(R, np.sqrt(f_nS))
		np.fill_diagonal(T, np.sqrt(f_Sp))
		np.fill_diagonal(D, np.sqrt(y))
		np.fill_diagonal(D_adj, np.sqrt(y_adj))

		# psi 		= np.dot(np.dot(np.dot(T , np.dot(Z.T, V)), D.I), D_adj)
		psi 		= T * (Z.I * V) * D.I * D_adj
		self.psi 	= psi.A

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
		# if len(phi.shape) == 1: #17/6/2016, removed as should be unneccesary if binarizing labels
		# 	C = 1
		# if len(phi.shape) > 1:
		# 	C = phi.shape[1]
		C = phi.shape[1]

		f_nS 	= self.f_nS
		f_Sp 	= self.f_Sp
		y 		= self.y
		y_adj 	= self.y_adj

		R 		= np.zeros((N,N))
		H 		= np.zeros((C,C))
		D 		= np.zeros((J,J))
		D_adj 	= np.zeros((J,J))

		R 		= np.mat(R)
		H 		= np.mat(H)
		D 		= np.mat(D)
		D_adj 	= np.mat(D_adj)
		phi		= np.mat(phi)
		V 		= np.mat(V)

		np.fill_diagonal(R, np.sqrt(f_nS))
		np.fill_diagonal(H, phi.sum(axis=0))
		np.fill_diagonal(D, np.sqrt(y))
		np.fill_diagonal(D_adj, np.sqrt(y_adj))

		# omega = np.dot(np.dot(np.dot(np.dot(np.dot(R, phi), H).T, V), D.I), D_adj)
		omega = ((R * phi * H).T * V) * D.I * D_adj
		self.omega = omega.A

	def calc(self):
		"""Summary
		Carries out mca analysis. User only need interact with this class method after instantiation of a mca class.

		"""
		self.decompose()
		if self.row:
			self.proj_row()
		if self.col:
			self.proj_col()
		if self.phi is not None:
			self.proj_feat()

	def evec(self, mode=None):
		if mode is None:
			return np.fliplr(self.V)
		else:
			return np.fliplr(self.V)[:,mode]
	def eval(self, mode=None):
		if mode is None:
			return np.flipud(self.y_adj)
		else:
			return np.flipud(self.y_adj)[mode]
	def eval_raw(self, mode=None):
		if mode is None:
			return np.flipud(self.y)
		else:
			return np.flipud(self.y)[mode]
	def score_row(self, mode=None):
		if mode is None:
			return np.fliplr(self.theta)
		else:
			return np.fliplr(self.theta)[:,mode]
	def score_col(self, mode=None):
		if mode is None:
			return np.fliplr(self.psi)
		else:
			return np.fliplr(self.psi)[:,mode]
	def score_feat(self, mode=None):
		if mode is None:
			return np.fliplr(self.omega)
		else:
			return np.fliplr(self.omega)[:,mode]

	def findI(self, alpha=0.01, use_continuity=True):#, zero_method='wilcox', correction=False):
		"""Summary
		Use the Wilcoxon rank-sum test to find the first I principle axes that contribute significant information.
		The method will properly re-order principle axes assuming the variables were not modified by the end-user.

		Parameters
		----------
		alpha
			critical p-value for determining if collection of columns is significant
		use_continuity
			parameter for scipy.stats.mannwhitneyu. if true, the method performs a continuity correction of 1/2.

		"""
		J = self.J
		N = self.N
		theta = np.fliplr(self.theta)
		y = np.flipud(self.y_adj)

		# Test to set max I is not working properly
		# lhs = np.empty((J,))
		# rhs = np.empty((J,))
		# for j in xrange(J):
		#     if j != 0:
		#         lhs[j] = y[0:j].sum()
		#         rhs[j] = (theta**2).sum(axis=0)[0:j].sum()
		#     if j == 0:
		#         lhs[j] = y[j]
		#         rhs[j] = (theta ** 2).sum(axis=0)[j]
		# diff = np.abs(lhs-rhs)
		# Imax = np.max([np.argmin(diff)+1, 2])
		# Imax = np.min([Imax, len(lhs)])

		pval = np.ones((len(range(1,J)),))
		for i in xrange(1,J):
			x = (theta**2)[:,0:i-1].sum(axis=1)
			y = (theta**2)[:,0:i].sum(axis=1)
			# results = stats.ranksums(x,y) # Other similar methods that are not preferred
			# results = stats.wilcoxon(x, y, zero_method=zero_method, correction=correction)
			# pval[i-1] = stats.norm.cdf(results[0])
			results = stats.mannwhitneyu(x, y, use_continuity=use_continuity)
			pval[i-1] = results[1]

		mask = (pval < alpha)
		I = np.where(mask==False)[0][0]
		self.I = I
		self.pval = pval
		return I

	def cluster(self, k, n_iter=10):
		I = self.I
		theta = np.fliplr(self.theta)[:,0:I]

		centroid, group = kmeans2(theta, k = k, iter=n_iter, thresh=1e-05, minit='random', missing='warn', check_finite=True)

		self.group = group
		self.centroid = centroid


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
			# axis = np.argmax(v.shape)
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