'''
Created on May 11, 2009

@author: lunt
'''

def flatten(x,maxDepth=-1):
	"""BLATENTLY STOLEN FROM : http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
	flatten(sequence) -> list
	Returns a single, flat list which contains all elements retrieved
	from the sequence and all recursively contained sub-sequences
	(iterables).

	Examples:
	>>> [1, 2, [3,4], (5,6)]
	[1, 2, [3, 4], (5, 6)]
	>>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
	[1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]
	"""
	result = []
	for el in x:
		#if isinstance(el, (list, tuple)):
		if (not maxDepth==0) and hasattr(el, "__iter__") and not isinstance(el, basestring):
			result.extend(flatten(el,maxDepth=(maxDepth-1)))
		else:
			result.append(el)
	return result

