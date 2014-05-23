r"""
Link class
"""
#*****************************************************************************
#  Copyright (C) 2014
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.free_group import FreeGroupElement
from sage.groups.braid import Braid
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

class Link:
"""
The base class for Link, taking input in three formats namely Briadword, gauss_code, dt_code
"""
	def __init__(self, input = None, gauss_code = None, dt_code = None):
		if type(input) == Braid:
			self._braid = input
			self._gauss_code = None
			self._dt_code = None

		elif gauss_code != None:
			self._braid = None
			self._gauss_code = gauss_code
			self._dt_code = None

		elif dt_code != None:
			self._braid = None
			self._gauss_code = None
			self._dt_code = dt_code

		else:
			raise Exception("Invalid input")

	def braidword(self):
		r"""
		Returns the braidword

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Braidword representation of the INPUT

		EXAMPLES::
		"""

		if self._braid != None:
			return list(self._braid.Tietze())

		elif self._gauss_code != None:
			return "Not implemented Error"

		elif self._dt_code != None:
			return "Not Implemented Error"

	def braid(self):
		r"""
		Returns the braid

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Braid representation of the INPUT

		EXAMPLES::
		"""

		if self._braid != None:
			return self._braid

		elif self._gauss_code != None:
			return "Not implemented Error"

		elif self._dt_code != None:
			return "Not Implemented Error"

	def gauss_code(self):
		r"""
		Returns the gauss_code

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Gauss code representation of the INPUT

		EXAMPLES::
		"""
		if self._gauss_code != None:
			return self._gauss_code

		elif self._braid != None:
			return "Not Implemented Error"

		elif self._dt_code != None:
			return "Not Implemented Error"

	def dt_code(self):
		r"""
		Returns the dt_code

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- DT Code representation of the INPUT

		EXAMPLES::
		"""
		if self._dt_code != None:
			return self._dt_code

		elif self._braid != None:
			return "Not Implemented Error"

		elif self._gauss_code != None:
			return "Not Implemented Error"

	def braidcomponentsmatrix(self):
		r"""
		Returns the braid components matrix # fails for input like 1 4 1 4 1 5 1 6

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Matrix containing the components is returned

		EXAMPLES::

			sage: from sage.knots import link
			sage: B = BraidGroup(4)
			sage: L = link.Link(B([-1, 3, 1, 3]))
			sage: L.braidcomponentsmatrix()
			[-1  1  0  0]
			[ 3  3  0  0]

		"""
		b = self.braid()
		ml = list(b.Tietze())
		l = [abs(k) for k in ml]
		sorting1 = list(set(range(l[0],max(l)+1)) - set(l))
		sorting = list(set(range(l[0],max(l)+1)) - set(sorting1))
		missing = []
		difference = []
		difference.append(sorting[0] -1)
		for i in range(1, len(sorting)):
			difference.append(sorting[i] - sorting[i-1])
		missing = list(set(range(sorting[0],sorting[-1]+1)) - set(sorting))
		if len(missing) == 0:
			A = matrix(ZZ,1,len(ml),ml)
			return A
		else:
			A = matrix(ZZ, len(missing) + 1, len(l))
			for i in range(len(missing) + 1):
	        		k = 0
	        		for j in range(len(l)):
	        			if i == 0:
	        				if (abs(l[j]) < missing[i]):
	        					A[0, k] = ml[j]
	       						k = k + 1
	       				elif i>len(missing):
	       					if missing[i] < abs(x[j]):
	        					A[i, k] = ml[j]
	        					k = k + 1
	        			else:
	        				if(missing[i-1] < abs(l[j])):
	       						A[i, k] = ml[j]
	 	       					k = k + 1
	 	
	 		for i in range(len(missing) + 1):
					l1 = list(A.row(i))
	 				M = min(l1)
	 				if M == 0:
	 					a = min(n for n in l1 if n!=M)
	 					for j in range(len(l)):
	 						if(A[i,j] != a and A[i,j] != a+1):
								A[i,j] = 0

			value = None
			for i in range(len(missing) + 1):
					for k in range(i+1, len(missing) + 1):
							if(A[i,0] == A[k,0]):
								value = 1
								A1 = A.delete_rows([k])
					return A1 if value == 1 else A				
			
	def braidwordcomponentsvector(self):
		r"""
		From braidcomponentsmatrix it is converted to a vector

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Vector containing non-zero values

		EXAMPLES::

			sage: from sage.knots import link
			sage: B = BraidGroup(4)
			sage: L = link.Link(B([-1, 3, 1, 3]))
			sage: L.braidwordcomponentsvector()
			[-1, 1, 3, 3]
		"""
		A = self.braidcomponentsmatrix()
		if(A.nrows()>0):
			x2 = []
			for i in range(A.nrows()):
				for j in range(A.ncols()):
					x2.append(A[i,j])
					x3 = [value for value in x2 if value != 0]
			return x3
		else:
			x3 = A.row(1)
			return list(x3)
	
	def homology_generators(self):
		r"""
		Returns the homology generators

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- The homology generators relating to the braid word representation

		EXAMPLES::

			sage: from sage.knots import link
			sage: B = BraidGroup(4)
			sage: L = link.Link(B([-1, 3, 1, 3]))
			sage: L.homology_generators()
			[1, 0, 3]
		"""
		x4 = self.braidwordcomponentsvector()
		hom_gen = []
		for j in range(len(x4)-1):
			a = abs(x4[j])
			for i in range(j+1, len(x4)):
					if(a == abs(x4[i])):
						hom_gen.append(i)
						break
			else:
				hom_gen.append(0)
		return hom_gen

	def Seifert_Matrix(self):
		r"""
		Returns the Seifert Matrix

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Returns the Seifert Matrix of the link.

		EXAMPLES::

			sage: from sage.knots import link
			sage: B = BraidGroup(4)
			sage: L = link.Link(B([-1, 3, 1, 3]))
			sage: L.Seifert_Matrix()
			[ 0  0]
			[ 0 -1]
		"""
		x5 = self.braidwordcomponentsvector()
		h = self.homology_generators()
		hl = len(h)
		A = matrix(ZZ, hl, hl)
		for i in range(hl):
			if h[i] != 0:
				for j in range(i,hl):
						if i == j:
							A[i,j] = -cmp((x5[i] + x5[h[i]]),0)
						elif (h[i] > h[j]):
							A[i,j] = 0
							A[j,i] = 0
						elif (h[i] <  j):
							A[i,j] = 0
							A[j,i] = 0
						elif (h[i] == j):
							if(x5[j] > 0):
								A[i,j] = 0
								A[j,i] = 1
							else:
								A[i,j] = -1
								A[j,i] = 0
						elif abs(abs(x5[i]) - abs(x5[j])) > 1:
							A[i,j] =  0
						elif (abs(x5[i]) - abs(x5[j]) == 1):
							A[i,j] = 0
							A[j,i] = -1
						elif (abs(x5[j])- abs(x5[i]) == 1):
							A[i,j] = 1
							A[j,i] = 0
						else:
							A[i,j] = 2
							A[j,i] = 2
			else:
				for k in range(hl):
					A[k,i] = 0
					A[i,k] = 0
		k = []
		for i in range(hl):
				if h[i] == 0:
					k.append(i)
		for i in reversed(k):
				A = A.delete_rows([i])
				A = A.delete_columns([i])
		return A


	def link_number(self):
		r"""
		Returns the link number

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Link number of the link

		EXAMPLES::

			sage: from sage.knots import link
			sage: B = BraidGroup(4)
			sage: L = link.Link(B([-1, 3, 1, 3]))
			sage: L.link_number()
			4
		"""
		count = 0
		b = self.braid().permutation()
		c = [0 for i in range(len(b))]
		for i in range(len(b)):
			if c[i] == 0:
				c[i] = 1
				k = b[i]
				while((k-1) != i):
					c[(k-1)] = 1
					k = b[k-1]
				count = count + 1
		return count


	def smallest_equivalent(self):
		r"""
		Returns the braidword

		INPUT:
			- Either a briadword, gauss_code, dt_code

		OUTPUT:
			- Smallest equivalent of the given braid word representation.

		EXAMPLES::

			sage: from sage.knots import link
			sage: B = BraidGroup(4)
			sage: L = link.Link(B([-1, 3, 1, 3]))
			sage: L.smallest_equivalent()
			[-1, 3, 1, 3]
		"""
		b = list(self.braid().Tietze())
		b1 = min([abs(k) for k in b])
		for i in range(len(b)):
			if b[i] > 0:
				b[i] = b[i] - b1 + 1
			else:
				b[i] = b[i] + b1 - 1
		return b
