import logging
import sys
import os
import logging
import copy
import itertools
import networkx

from collections import defaultdict
import numpy as np
from collections import namedtuple

MAXl = 12

class Clique(object):
	def __init__(self, node_list, cid):
		self.ID = str(cid)
		self.node_list = sorted(node_list)
		self.node_set = set(node_list)
		self.clique_size = len(node_list)
		return
	
	#merge two clique
	def merge(self, otherC):
		self.node_set = self.node_set.union(otherC.node_set)
		self.node_list = list(self.node_list)
		self.ID = self.ID+otherC.ID
	#descending order
	def __cmp__(self, other):
		#return cmp(other.cliqueSize, self.cliqueSize)
		return cmp(self.clique_size, other.clique_size)
		
	def __str__(self):
		return str(self.ID) + ', '+ str(self.node_list) + ', with size: ' + str(len(self.node_list))

class CliqueList(object):
	def __init__(self, items = []):
		self.list = []
		clique_id = 1
		for elem in items:
			self.list.append(Clique(elem, clique_id))
			clique_id = clique_id+1

		self.list = sorted(self.list, key = lambda clique: clique.clique_size, reverse = True)			

	def add(self, item):
		self.list.append(item)
	
	def toPureList(self):
		result = []
		for elem in self.list:
			result.append(elem.node_list)
		
		return result

	def merge_close_cliques(self):
	
		maxScore = -sys.maxint - 1
		bestk2 = -1
		bestj = -1
		clique_dict={}
		
		pass_list = range(0, len(self.list))
			
		result_dict = {}
		
		i = pass_list[0]
		while i in pass_list:
			
			k1 = self.list[i].ID
			c1 = self.list[i].node_set
			print self.list[i]
			for j in pass_list:
				k2 = self.list[j].ID
				c2 = self.list[j].node_set
				
				if k1 != k2:
					score = len(c1.intersection(c2))
					if score > maxScore:
						maxScore = score
						bestk2 = k2
						bestj = j

			if bestj != -1:
				c1 = c1.union(c2)
				key = k1 + bestk2
				result_dict[key] = Clique(list(c1), key)
				pass_list.remove(i)
				pass_list.remove(bestj)
				
				bestk2 = -1
				bestj = -1
				maxScore = -sys.maxint - 1

			if len(pass_list) <= 1:
				break
				
			i = pass_list[0]
					
		return result_dict

	
	def __str__(self):
		outputStr = ''
		for elem in self.list:
			outputStr = outputStr + str(elem) + '\n'
		return outputStr

class Basis(object):
	def __init__(self):
		self.ID = ''
		self.node_list = []
		self.node_set = set()
		self.basis_size = 0
		return
	
	def update2(self, nodeList):
		self.ID = '0'
		self.node_list = nodeList
		self.node_set = set(nodeList)
		self.basis_size = len(nodeList)
		return

	def update(self, clique):
		self.ID = clique.ID # use the clique ID to form the basis ID, both of them are strings
		self.node_list = clique.node_list
		self.node_set = set(clique.node_list)
		self.basis_size = len(self.node_list)
		return
	
	#merge the clique to this basis
	def merge(self, clique):
		self.ID = self.ID + clique.ID
		self.node_set = self.node_set | clique.node_set
		self.node_list = list(self.node_set)
		self.basis_size = len(self.node_list)

	#merge the clique to this basis
	def merge_node_list(self, node_list):
		self.node_set = self.node_set | set(node_list)
		self.node_list = list(self.node_set)
		self.basis_size = len(self.node_list)

		
	def __str__(self):
		return 'set(' + str(self.node_list) + '), '
		
	def calcIntersection(self, otherBasis):
		return len(self.node_set.intersection(otherBasis.node_set))

	def calcUnion(self, otherBasis):
		return len(self.node_set.union(otherBasis.node_set))

# The container of a set of basis. This also forms a solution. 
class Bases(object):
	def __init__(self):
		self.basis_dict = {}
		self.l = 0
		self.w = 0
		self.score = 0
		self.ID = ''
		self.candi_list = []
		self.optL = 6
		self.optL_2 = 3
		self.maxID = 0
		self.top_eta = []
		self.top_eta_pairs = []
		self.merge_count = 0
		self.prev_merge_count = -1		
		
	def toList(self):
		basis_list = []
 		for k, b in sorted(self.basis_dict.items()):
			basis_list.append(b)
			
		return basis_list
			
	def addSingleItems(self, orphans):

		start = 0
		for k, b in sorted(self.basis_dict.items()):
			if b.basis_size < self.optL:
				length = self.optL - b.basis_size
				self.basis_dict[k].merge_node_list(orphans[start: start + length])
				start = start + length
		
		i = 1
		item_list = []
		count = 1 
		
		for item in orphans[start:]:
			item_list.append(item)
			if i % self.optL_2 == 0:
				newBasis = Basis()
				newBasis.update2(item_list)
				self.basis_dict[str(int(self.maxID) + count)] = newBasis
				count = count + 1
				item_list = []
		
			i = i + 1
		
		if len(item_list) == 0:
			return
		else:
			newBasis = Basis()
			newBasis.update2(item_list)
			self.basis_dict[str(int(self.maxID) + count)] = newBasis
		
		
	def checkStopping(self, optL):
		if self.prev_merge_count < self.merge_count:
			self.prev_merge_count = self.merge_count
			return True
		
		if self.prev_merge_count == self.merge_count:
			return False
			
	def make_orphan_node_list(self, single_list, clique_list):
		
		single_set = set(single_list)
		clique_node_list = make_node_list_2(clique_list)
		clique_node_set = set(clique_node_list)
		orphan_list = list(single_set.difference(clique_node_set))
		return orphan_list
	
	#the main logic of algo2 
	def optimalMerge(self, atop_eta, atop_eta_pairs, cliqueList, aoptL):
		self.optL = aoptL
		prevScore = sys.maxint
		stopKeys = []
		
		self.top_eta = atop_eta
		self.top_eta_pairs = atop_eta_pairs
		
		while self.checkStopping(self.optL):

			maxIntersection = -sys.maxint - 1
			bestk1 = -1
			bestk2 = -1
			finalUnionScore = 0
			
			prevScore = self.score		
			for k1, b1 in sorted(self.basis_dict.items()):
				if k1 in stopKeys:
					continue
				if b1.basis_size >= self.optL:
					stopKeys.append(k1)
					continue
				for k2, b2 in sorted(self.basis_dict.items()):
					if k2 in stopKeys:
						continue

					if k2 != k1:
						intersectionScore = b1.calcIntersection(b2)
						unionScore = b1.calcUnion(b2)
						if intersectionScore > maxIntersection:
							if unionScore <= self.optL:
								maxIntersection = intersectionScore
								finalUnionScore = unionScore
								bestk1 = k1
								bestk2 = k2
							else:
								maxIntersection = -sys.maxint - 1
								bestk1 = -1
								bestk2 = -1
								finalUnionScore = 0
								continue
			
			if bestk1 != -1 and bestk2 != -1:
				self.merge2(bestk1, bestk2)
				
		clique_list_list = cliqueList.toPureList()
		orphans = self.make_orphan_node_list(self.top_eta, clique_list_list)
		
		self.addSingleItems(sorted(orphans))
		self.update_score()
		
	# calculate how many items has been covered by this set of basis
	def calc_coverage(self):
		items_covered = set()
		for key, basis in self.basis_dict.iteritems():
			items_covered = items_covered.union(set(basis.node_list))
		
		top_eta_set = set(self.top_eta)
		return len(items_covered)

	#make the clique as a new basis and add it to the basis list
	def appendClique(self, clique):
		newBasis = Basis()
		newBasis.update(clique)
		self.basis_dict[clique.ID] = newBasis
		self.update_score()
	
	#make the clique as a new basis and add it to the basis list
	def append(self, clique, targetID):
		if targetID == '0':
			newBasis = Basis()
			newBasis.update(clique)
			self.basis_dict[clique.ID] = newBasis
			del self.basis_dict['0']
		else:
			targetBasis = self.basis_dict[str(targetID)]
			oldBasisID = targetBasis.ID
			targetBasis.merge(clique)
			self.basis_dict[targetBasis.ID] = targetBasis
			del self.basis_dict[oldBasisID]
		self.update_score()
	
	# add an existing basis to the basis list
	def unite(self, basis):
		self.basis_dict[basis.ID] = basis
		self.update_score()
	
	def update_score(self):
		lenList = [value.basis_size for key, value in self.basis_dict.iteritems()]	
		self.l = max(lenList)
		self.w = len(self.basis_dict)
		self.score = self.w**2 * 2**self.l
		self.ID = ''
		for key, value in self.basis_dict.iteritems():
			self.ID = self.ID + value.ID
		
		intKeys = [int(elem) for elem in self.basis_dict.keys()]
		self.maxID = max(intKeys)
	
	def merge(self, clique):
		# merge it to all the existing ones
		for basis in self.basis_list:
			basis.merge(clique)
		self.update_score()

	#merge with another basis in the basis_dict
	def merge2(self, k1, k2):
		self.merge_count = self.merge_count + 1
		b1 = self.basis_dict.get(k1)
		b2 = self.basis_dict.pop(k2)
		b1.merge(b2)
		self.update_score()
		

	def __str__(self):
		outputStr = ''
		outputStr = outputStr + '\n(l = ' + str(self.l) + ', w = ' + str(self.w) + '), score = ' + str(self.score) + ', coverage = ' + str(self.calc_coverage()) + '\n'
		for key, basis in self.basis_dict.iteritems():
			outputStr = outputStr + '\t' + str(basis) + '\n'
		return outputStr

	#descending order
	def __cmp__(self, other):
		#return cmp(other.cliqueSize, self.cliqueSize)
		return cmp(self.score, other.score)



def networkx_process(top_node_list, top_pair_list):
		
	g = networkx.Graph()
	g.add_nodes_from(top_node_list)
	g.add_edges_from(top_pair_list)
	clique_list = list(networkx.find_cliques(g))
	
	return clique_list

# make the nodes involved in the pairs
def make_node_list(pair_list):
	top_node_set = set()
	for pair in pair_list:
		top_node_set.add(pair[0])
		top_node_set.add(pair[1])
	top_node_list = []
	top_node_list = list(top_node_set)
	return top_node_list
# [[], []]
# make the nodes involved in the cliques
def make_node_list_2(clique_list):
	node_set = set()
	for clique in clique_list:
		for elem in clique:
			node_set.add(elem)
	return list(node_set)

def get_optimal_bases(single, pairs):	
	top_single_list = single
	top_pair_list = pairs
	top_node_list = make_node_list(top_pair_list)
	clique_list = networkx_process(top_node_list, top_pair_list)
	
	clique_list.sort(key=len)
	
	optL = len(clique_list[-1])
	
	cliqueList = CliqueList(clique_list)
	bases = Bases()
	for clique in cliqueList.list:
		bases.appendClique(clique)

	bases.optimalMerge(single, pairs, cliqueList, optL)
	
	return bases
	
	
def main():
	logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
	logging.basicConfig(level=logging.INFO, stream=sys.stdout)

	from time import clock, time

	starttime = time()

	bases = get_optimal_bases(top_eta, top_eta_pairs)
	print 'bases = ', bases

	endtime = time()

	interval= endtime - starttime
	print 'time = ', str(interval)

if __name__ == '__main__':
	main()



