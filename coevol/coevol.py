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
	from Bio import SeqIO, Seq, Alphabet
	from scipy import linalg, stats
	from scipy.cluster.vq import kmeans2
except:
	print('ERROR: required libraries could not be loaded')


"""
define coevol classes
"""
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
		row :
			project rows onto principle axes if True
		col :
			project columns onto principle axes if True
		tol:
			eigenvalues falling within the range of 1-tol<y<1+tol will be discarded as trivial

		"""
		self.phi, self.key_lbl = binarize(phi, axis=1, lbl=True)
		self.row = row
		self.col = col
		self.tol = np.abs(tol)

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
		"""
		Summary
		Quickly return first eigenvector according to the mode selected
		"""
		if mode is None:
			return np.fliplr(self.V)
		else:
			return np.fliplr(self.V)[:,mode]
	def eval(self, mode=None):
		"""
		Summary
		Quickly return first eigenvalue according to the mode selected
		"""
		if mode is None:
			return np.flipud(self.y_adj)
		else:
			return np.flipud(self.y_adj)[mode]
	def eval_raw(self, mode=None):
		"""
		Summary
		Quickly return first corrected eigenvalue according to the mode selected
		"""
		if mode is None:
			return np.flipud(self.y)
		else:
			return np.flipud(self.y)[mode]
	def score_row(self, mode=None):
		"""	
		Summary
		Quickly return the row projection corresponding to the selected mode
		"""
		if mode is None:
			return np.fliplr(self.theta)
		else:
			return np.fliplr(self.theta)[:,mode]
	def score_col(self, mode=None):
		"""	
		Summary
		Quickly return the column projection corresponding to the selected mode
		"""
		if mode is None:
			return np.fliplr(self.psi)
		else:
			return np.fliplr(self.psi)[:,mode]
	def score_feat(self, mode=None):
		"""	
		Summary
		Quickly return the feature projection corresponding to the selected mode
		"""
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
		I = np.where(mask==False)[0][0] # First occurance of False indicates number of principle axes to be used!
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

def filter_col_gaps(data, cutoff=0.1):
	"""Summary
	Function accepts a sequence alignment and removes all columns (sequence positions) where gaps occur >= [cutoff] frequency.

	Parameters
	----------
	data : ndarray
		ndarray from numpy containing sequence alignment (columns are positions, rows are sequences)
	cutoff : float
		if gaps occur greater than or equal to this cutoff frequence, then the column will be removed
	"""
	mask = (data == '-')
	n_row = float(data.shape[0])
	count = mask.astype(int).sum(axis=0)
	freq = count / n_row
	mask_keep = freq < cutoff
	return data[:, mask_keep]


def parse_aln(filein, style='fasta'):
	"""Summary
	Function parses a sequence alignment into a numpy array

	Parameters
	----------
	filein : str
		full path to file containing sequence alignment
	style : str
		format of the sequence alignment (only 'fasta' has been tested, other formats compatible with Biopython should work)
	"""
	obj = SeqIO.parse(filein, style)
	fa = [x for x in obj]
	n_pos = np.max([len(x.seq) for x in fa])
	n_seq = len(fa)
	arr = np.empty((n_seq, n_pos)).astype('str')
	for i in xrange(n_seq):
		arr[i,] = list(fa[i].seq)
	return arr


def write_aln(ids, arr, fileout, style='fasta'):
	"""Summary
	Function will write a sequence alignment to file from numpy arrays containing the ids and sequences

	Parameters
	----------
	ids : list or ndarray
		Array of ids corresponding to each row of sequence array. If provided as a list, it will be coerced to an ndarray.
	arr : list or ndarray
		Array of sequences (one character per element). If provided as a list, it will be coerced to an ndarray.
	fileout : str
		Full path where sequence alignment should be written
	style : str
		String specifying format of style to be written. Must agree with options from Biopython. (Currently only "fasta" tested)
	"""
	if isinstance(ids, np.ndarray) == False:
		ids = np.asarray(ids)
	if isinstance(arr, np.ndarray) == False:
		ids = np.asarray(ids)
	aln = []
	for x, y in zip(ids.tolist(), arr.tolist()):
		seq = ''.join(y)
		seq = Seq.Seq(seq, Alphabet.SingleLetterAlphabet())
		id = x
		aln.append(SeqIO.SeqRecord(seq=seq, id=id, name=id, description=''))
	SeqIO.write(aln, fileout, style)


def parse_lbl(filein, schar=None, index=0, style='fasta', sfmt='|S128'):
	"""Summary
	Function to parse a label designation from a header in a sequence alignment (must be in id field).

	Parameters
	----------
	filein : str
		file containing sequence alignment
	schar : str or None
		character to split header on. if None, headers are not split and entire header is returned.
	index : int or None
		python index for element of the list (following string splitting) containing desired label.
	style : str
		string for format of file to be read (must be compatible with Biopython)
	sfmt : str
		format for ndarray where labels will be placed. if size is too small or too large for entire label, increase by specifying size here
	"""
	obj = SeqIO.parse(filein, style)
	fa = [x for x in obj]
	n_seq = len(fa)
	arr = np.empty((n_seq,), dtype='|S128')
	for i in xrange(n_seq):
		if schar is None:
			arr[i] = fa[i].id
		else:
			arr[i] = fa[i].id.split(schar)[index]
	return arr
