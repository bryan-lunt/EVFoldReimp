'''
Make CNS constraints for contacts.


Created on Feb 12, 2013

@author: lunt
'''

from itertools import izip,count

constraintstr = 'assign (resid {0:n} and name {2:s})  (resid {1:n} and name {3:s})  {4:.1f} {5:.1f} {6:.1f} weight {7:.5f}\n'
EVFold_Distances = (4.0,4.0,3.0)

#See EVFold Appendix A10
res_specific_atom = {'A':'CB','C':'SG','D':'OD1',
					'E':'OE1','F':'CZ','G':'CA',
					'H':'CE1','I':'CD1','K':'NZ',
					'L':'CD1','M':'CE','N':'OD1',
					'P':'CG','Q':'OE1','R':'NH1',
					'S':'OG','T':'OG1','V':'CG1',
					'W':'CH2','Y':'OH'}


def make_contact_constraints(pairs_and_weights,one_sequence,distances=EVFold_Distances):
	
	returnStr = ''
	
	for i,j,weight in pairs_and_weights:
		#CA constraint
		returnStr += constraintstr.format(*((i+1,j+1,'CA','CA') + distances + (weight,)))
		#CB constraint
		if not(one_sequence[i] == 'G' or one_sequence[j] == 'G'):
			returnStr += constraintstr.format(*((i+1,j+1,'CB','CB') + distances + (weight,)))
		#residue specific
		i_atom = res_specific_atom[one_sequence[i]]
		j_atom = res_specific_atom[one_sequence[j]]
		returnStr += constraintstr.format(*((i+1,j+1,i_atom,j_atom) + distances + (weight,)))

	return returnStr


def make_EVFold_contacts(pairs,one_sequence):
	
	EVFold_startweight = 10.0
	
	weightedlist = [ (i,j,EVFold_startweight/num) for (i,j),num in izip(pairs,count(1))]

	return make_contact_constraints(weightedlist,one_sequence,distances=EVFold_Distances)
	