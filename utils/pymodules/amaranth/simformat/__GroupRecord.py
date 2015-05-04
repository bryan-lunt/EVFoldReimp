'''
Created on Sep 21, 2010

@author: lunt
'''
import Bio.SeqRecord

class MultiDictView(object):
	def __init__(self,dicts=[]):
		self.dicts = dicts
		
	def __getitem__(self,key):
		retval = [i.get(key,None) for i in self.dicts]
		if all([i is None for i in retval]):
			raise KeyError(key)
		return tuple(retval)
	
	def get(self,key,*varargs):
		default = None if len(varargs) == 0 else varargs[0]
		retval = [i.get(key,default) for i in self.dicts]
		return tuple(retval)
	
	def has_key(self,key):
		return any([i.has_key(key) for i in self.dicts])


class SeqRecordGroup(Bio.SeqRecord.SeqRecord):
	'''
	This stores a list of regular SeqRecords, and presents them as a single joined record.
	'''

	def __init__(self,params):
		'''
		Constructor
		'''
		self.__mylist = list()
	
	
	def getID(self):
		return ':'.join([i.id for i in self.__mylist])
	
	def setID(self,newval):
		raise Exception("You Should Not Do That.")
	
	id = property(getID,setID)