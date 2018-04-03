import sys
import os
import logging
from random import shuffle

from collections import defaultdict

import numpy
from collections import namedtuple

from fp_growth import find_frequent_itemsets

class Dataset(object):
	
	def __init__(self, path, tau=None, maxlen = float('inf')):
		self.path = path
		self.tau = tau
		
		self.Itemset = namedtuple('Itemset', ['elem', 'sup'])
		
		self.transactions = []
		self.items = defaultdict(set)
		
		self.topk = None
		
		self.fi = None

		self.item_info = {} #key=item, value=sup, minsup=0
		
		self._load(path, maxlen)
		

		
		self.item_sup ={}
		
		for item in self.items:
			self.item_sup[item] = len(self.items[item])		
	
	def mine_freq_item(self, minsup, maxlen):
#		logging.debug('Mining frequent itemsets')
		self.fi = self.get_freq_itemsets(minsup, maxlen)
		self.topk = [item for item in self.items]
		self.topk.sort(key=lambda item: len(self.items[item]), reverse=True)
		
#		logging.debug('END Mining frequent itemsets')
	
	def revert(self,):
		self.tau = None
		self.transactions = self.original_trans
	
	def limit_contribution(self, tau):
		if self.tau == tau:
			return
		
		self.tau = tau
		new_transactions = []
		self.original_items = self.items
		self.items = defaultdict(set)
		
		throw = 0
		for i, trans in enumerate(self.transactions):
			
			throw += max(0, len(trans) - tau)
			shuffle(trans)
			trans = trans[:tau]
			
			for item in trans:
				self.items[item].add(i)
			
			new_transactions.append(trans)
		
		self.original_trans = self.transactions
		self.transactions = new_transactions
		
		self.item_sup ={}
		for item in self.items:
			self.item_sup[item] = len(self.items[item])
		
		return throw
	
	def psuedo_limit_contribution(self, tau):
		itemsup = defaultdict(int)
		
		throw = 0
		for i, trans in enumerate(self.transactions):
			throw += max(0, len(trans) - tau)
			
			f = min(tau, len(trans))*1.0/len(trans)#tau*1.0/len(trans)#min(tau, len(trans))*1.0/len(trans)

			for item in trans:
				itemsup[item] += f

		return throw, itemsup
		
	def psuedo_limit_contribution2(self, tau):
		itemsup = defaultdict(int)
		
		throw = 0
		for i, trans in enumerate(self.transactions):
			throw += max(0, len(trans) - tau)
			
			shuffle(trans)
			trans = trans[:tau]
			
			f = min(len(trans)*1.0/tau, 1.0)

			for item in trans:
				itemsup[item] += f

		return throw, itemsup
		
	def get_sup_of_kth_item(self, k):
		return len(self.items[self.topk[k]])
	
	def get_sup_of_kth_itemset(self, k):
#		print "k = ", k
#		print "len(self.fi) = ", len(self.fi)
#		print "m = ", self.m
#		print "self.fi = ", self.fi

		return self.fi[k].sup
	
	def get_actual_topk_singleton(self, k):
		return self.topk[:k]

	def get_actual_topk_mixed(self, k):
		return self.fi[:k]	
		
	def _load(self, path, maxlen):
		f = open(path, 'r')
		
		len_t = []
		for i, line in enumerate(f):
			t = []
			for item in line.split():
				item = int(item)
				t.append(item)
			
			if len(t) > maxlen:
				continue
				
			len_t.append(len(t))
			if self.tau != None:
				shuffle(t)
				t = t[:self.tau]
			
			for item in t:
				self.items[item].add(i)
				if item not in self.item_info:
					self.item_info[item] = 1
				else:
					self.item_info[item] += 1
			
			self.transactions.append(t)
		
		f.close()
		
		self.avg_t_len = numpy.mean(len_t)
		self.median_t_len = numpy.median(len_t)
		self.max_t_len = numpy.max(len_t)
		self.m = len(self.items)
		self.n = len(self.transactions)
		

	def exclude(self, itemsets):
		for itemset in itemsets:
			print itemset
			self.fi.remove(itemset)
	
	def get_ith_most_freq(self, i):
		try:
			return self.fi[i]
		except:
			print i, len(self.fi)
			raise
	
	
	

	def get_freq_itemsets(self, min_sup, max_len):
		
		fi = []
		if self.path == '../datasets-space/kosarak-full-space.data':
			f = open('kosarak-minsupp-0.6-percent.txt')
			
			for line in f:
				tokens = line.split(',')
				itemsetstr = tokens[0]
				itemsetstr = itemsetstr.strip('[]')
				
				itemset = []
				for item in itemsetstr.split():
					itemset.append(int(item))
				
				sup = int(tokens[1])
				if sup < min_sup:
					break
				
				if len(itemset) <= max_len:	
					fi.append(self.Itemset(itemset, sup))
			
			f.close()
			return fi


		if self.path == '../datasets-space/aol-full-space.data':
			f = open('aol-minsupp-0.2-percent.txt')
			
			for line in f:
			
				tokens = line.split()
				itemset = []
				for item in tokens[:-2]:
					itemset.append(int(item))
				
				sup = int(tokens[-2])
				
				if sup < min_sup:
					break
				
				if len(itemset) <= max_len:	
					fi.append(self.Itemset(itemset, sup))

			#print "test loading aol, fi = ", fi
						
			f.close()
			return fi
		
		if max_len == 1:
			x = 0
			for item, coverage in self.items.iteritems():
				sup = len(coverage)
				
				if sup >= min_sup:
					x+=1
					fi.append(self.Itemset([item,], sup))
		else:
			x = 0
			for itemset, sup in find_frequent_itemsets(self.transactions, min_sup, include_support=True):
				if len(itemset) > max_len:
					continue
				fi.append(self.Itemset(itemset, sup))
				x+=1
#				if x % 100 == 0:
					#logging.debug('%d frequent itemsets obtained'%x)
			
#		logging.debug('Total for min sup %d = %d frequent itemsets'%(min_sup,x))
		
		fi.sort(key = lambda(i):i.sup, reverse = True)
						
		return fi
	
	def print_summary(self):
		print 'Average Transaction Length', self.avg_t_len
		print 'Median Transaction Length', self.median_t_len
		print 'Max Transaction Length', self.max_t_len
		print 'm', self.m
		print 'n', self.n
		

def test_noise(b_list, n):
	
	for b in b_list:
		noise = numpy.random.laplace(loc=0.0, scale=b*1.0, size=n)
		noise = [abs(v) for v in noise]
		print b, numpy.mean(noise), numpy.median(noise)
		

def main():
	logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
	d = Dataset('../datasets-space/mushroom-anon-dataset/tildeD-1-0.1.txt')
	d.mine_freq_item(4865, 200)
	
if __name__ == '__main__':
	main()

