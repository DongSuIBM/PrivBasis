import logging
import sys
import os
import logging
import copy
import itertools
import networkx
import math
from collections import defaultdict
import numpy as np
from collections import namedtuple
from dataset import Dataset
from algo2 import get_optimal_bases, Basis
from scipy.stats.stats import sem as stderr

ALPHA_1 = 0.1
ALPHA_2 = 0.4
ALPHA_3 = 0.5
ETA = 1.1



class Algo3(object):

	def __init__(self, dataset, ak, aepsilon):
		super(Algo3, self).__init__()
		self.dataset = dataset
		self.k = ak
		self.epsilon = aepsilon
		self.bases = []
		self.item_info_sorted = []
		
	#sample j from {1,2,...,k}
	def getLambda(self):
		N = len(self.dataset.transactions)		
		k1 = int(math.ceil(self.k*ETA))
		theta = self.dataset.get_sup_of_kth_itemset(k1)
		lenI = len(self.dataset.items)
		p = []
		self.item_info_sorted = sorted(self.dataset.item_info.items(), lambda x, y: cmp(y[1], x[1]))

		for i in range(0, lenI):
			f = self.item_info_sorted[i][1]
			exponent = self.epsilon*0.5*(1 - abs(f - theta))
			p.append(math.e**exponent)
		
		sumP = sum(p)
		p = [elem*1.0/sumP for elem in p]
		r = np.random.random()
		tmpSum = 0.0
		result = -1
		for i in range(0, lenI):
			tmpSum += p[i];
			if tmpSum >= r: 
				result = i
				break

		return result + 1
	
	def privBasisMain(self):
	
		self.lambdaVal = self.getLambda()
	
		if self.lambdaVal <= 12:
			F = self.GetFreqItems(self.lambdaVal, ALPHA_2*self.epsilon)
			results = self.BasisFreq( (1 - ALPHA_1 - ALPHA_2)*self.epsilon)			
			return results
			
		else:
			lambda_2 = ETA*self.k - self.lambdaVal
			lambda_2 = lambda_2/math.sqrt(max(1, lambda_2/self.lambdaVal))
			beta_1 = ALPHA_2*self.lambdaVal/(self.lambdaVal + lambda_2)
			beta_2 = ALPHA_2 - beta_1
						
			F = self.GetFreqItems(self.lambdaVal, beta_1*self.epsilon)
			pairs = self.getTopPairs(F)
			P = self.GetFreqElements(pairs, lambda_2, beta_2*self.epsilon)
			#algo 2
			bases_obj = get_optimal_bases(F, P)
			self.bases = bases_obj.toList()
			results = self.BasisFreq( (1 - ALPHA_1 - ALPHA_2)*self.epsilon)
			return results

	def GetFreqElements(self, pairsDict, lambdaValue, epsilonValue):
		N = len(self.dataset.transactions)
		p = []
		U_copy = copy.deepcopy(pairsDict.keys())
		exponents_list = []
		
		for i in range(len(U_copy)):
			f = pairsDict[U_copy[i]]
			expo = f*epsilonValue/lambdaValue
			exponents_list.append(expo)
		#deal with the overflow problem
		exponents_list = np.array(exponents_list, float).clip(min(exponents_list), 500).tolist()
		for expo in exponents_list:
			p.append(math.e ** (expo))
		
		sumP = sum(p)
		p = [elem*1.0/sumP for elem in p]
		X = []
		sampStruct = self.packItemAndProb(U_copy, p)
		count = 0
		while count < lambdaValue:
			r = np.random.random()
			tmpSum = 0.0
			sample = None
			sampleIndex = -1
			p = [k2 for k1, k2 in sampStruct]
			sumP = sum(p)
			p = [elem*1.0/sumP for elem in p]
			for j in range(len(p)):
				tmpSum += p[j];
				if tmpSum >= r: 
					sampleIndex = j
					sample = sampStruct[j][0]
					break
			del sampStruct[sampleIndex]
			count +=  1
			X.append(sample)
		
		return X
	
	def GetFreqItems(self, lambdaValue, epsilonValue):
		N = len(self.dataset.transactions)	
		p = []
		U_copy = [k1 for k1, k2 in self.item_info_sorted]
		exponents_list = []
		for i in range(len(self.item_info_sorted)):
			f = self.item_info_sorted[i][1]
			expo = f*epsilonValue/lambdaValue
			exponents_list.append(expo)

		#deal with the overflow problem
		exponents_list = np.array(exponents_list, float).clip(min(exponents_list), 500).tolist()
		for expo in exponents_list:
			p.append(math.e ** (expo))
		
		sumP = sum(p)
		p = [elem*1.0/sumP for elem in p]
		X = []
		sampStruct = self.packItemAndProb(U_copy, p)
		
		count = 0
		while count < lambdaValue:
			r = np.random.random()
			tmpSum = 0.0
			sample = None
			sampleIndex = -1
			
			p = [k2 for k1, k2 in sampStruct]
			sumP = sum(p)
			p = [elem*1.0/sumP for elem in p]

			for j in range(len(p)):
				tmpSum += p[j];
				if tmpSum >= r: 
					sampleIndex = j
					sample = sampStruct[j][0]
					break
			
			del sampStruct[sampleIndex]
			count +=  1
			X.append(sample)
				
		newBasis = Basis()
		newBasis.update2(X)
		self.bases.append(newBasis)
		
		return X
	
	def packItemAndProb(self, items, p):
		sampStruct = []
		for i in range(len(items)):
			sampStruct.append((items[i], p[i]))
		
		return sampStruct
	
	def getTopPairs(self, F):
		pairs_support = {}
		for i in xrange(self.lambdaVal):
			icoverage = self.dataset.items[F[i]]
			for j in xrange(i+1, self.lambdaVal):
				jcoverage = self.dataset.items[F[j]]
				pairs_support[(F[i], F[j])] = len(icoverage.intersection(jcoverage))

		return pairs_support

	def BasisFreq(self, eps):
		itemsets_all = defaultdict(list)
		w = len(self.bases)
		for basis in self.bases:
			l = len(basis.node_list)
			for itemset, (count, ncount) in self.reconstruct(basis.node_list, eps/w).iteritems():
				itemset = list(itemset)
				itemset.sort()
				itemset= tuple(itemset)
				itemsets_all[itemset].append((count, ncount, l))
		
		itemsets = [(itemset, countlist[0][0], self.weighted_avg(countlist)) for (itemset, countlist) in itemsets_all.iteritems()]
		itemsets.sort(key=lambda (i,c,nc):nc, reverse=True)
		
		return itemsets[:self.k]


	def reconstruct(self, basis, eps):
		
		bin_counts = {}
		
		for bin in self._get_subsets(basis):
			coverage = copy.copy(self.dataset.items[bin[0]])
			for item in bin[1:]:
				coverage.intersection_update(self.dataset.items[item])
			
			for item in basis:
				if item not in bin:
					coverage.difference_update(self.dataset.items[item])
			
			supp = len(coverage)
			bin_counts[tuple(bin)] = [supp, supp+np.random.laplace(scale=1.0/eps)]
		
		itemset_counts = {}
		for itemset in bin_counts:
			itemset_counts[itemset] = [0,0]
			for bin in bin_counts:
				if set(itemset).issubset(bin):
					bin = tuple(bin)
					itemset_counts[itemset][0] += bin_counts[bin][0]
					itemset_counts[itemset][1] += bin_counts[bin][1]
		
		return itemset_counts


	def _get_subsets(self, basis):
		n = len(basis)
		
		for index in xrange(1, 2**n):
			subset = []
			for i in xrange(n):
				if (index>>i)&1==1:	
					subset.append(basis[i])
					
			yield subset

	def weighted_avg(self, count_list):
		return np.mean([count_list[i][1] for i in xrange(len(count_list))])


class ExperimentAlgo3(object):
	
	def __init__(self, dataset):
		self.dataset = dataset
	
	def calc_accuracy_stat(self, actualset, expset):
		TP = len(actualset.intersection(expset))
		FN = len(actualset) - TP
		FP = len(expset.difference(actualset))

		fnr = FN*1.0/(len(actualset))
		precision = TP*1.0/(TP + FP)
		recall = TP*1.0/(TP+FN)

		return {'fnr': fnr, 'precision': precision, 'recall': recall}
	
	def run_experiment(self, eps, k, initial_minsup=None, repeat=5):
		# get accurate result
		if initial_minsup:
			self.dataset.mine_freq_item(initial_minsup, float('inf'))
		accurate_topk = self.dataset.fi[:k]
		
		#print accurate_topk
		accurate_itemsets = set([tuple(sorted(itemset.elem)) for itemset in accurate_topk])
		
		# get anonymized results
		results = defaultdict(list)
		for r in xrange(repeat):
			exp = Algo3(self.dataset, k, eps)
			anonymized_topk = exp.privBasisMain()
			#print anonymized_topk
		
			anonymized_itemsets = set(itemset for (itemset, c, nc) in anonymized_topk)
			rel_error = [abs(count - ncount)*1.0/count if count > 0 else 0 for (i, count, ncount) in anonymized_topk]
		
			stats = self.calc_accuracy_stat(accurate_itemsets, anonymized_itemsets)
			
			for st, val in stats.iteritems():
				results[st].append(val)
			results['relative_err'].append(np.median(rel_error))
		
		return results	
		
def main():
	logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
	logging.basicConfig(level=logging.INFO, stream=sys.stdout)
	
	datapath = {
		#'mushroom': ('../datasets-space/mushroom-full-space.data', 3300),
		'pumsbstar': ('../datasets-space/pumsb_star-full-space.data', 24000),
		# 'kosarak': ('../datasets-space/kosarak-full-space.data', 8000),
		# 'retail': ('../datasets-space/retail-full-space.data', 750),
		# 'AOL': ('../datasets-space/aol-full-space.data', 12000)				
	}

	#eps_list = [0.25, 0.5, 0.75, 1.0,]
	eps_list = [0.05, 0.1, 0.2, 0.25, 0.4, 0.5, 0.75, 0.8, 1.0]
	
	experiment = {
			
			#'mushroom': [50, 100], #e.g. 50 and 100 are k values
			#'mushroom': [100], #e.g. 50 and 100 are k values
			'pumsbstar': [100],
			#'kosarak': [100, 200, 300, 400],
			#'retail': [50, 100],
			#'AOL': [100, 200],			
	}
	
	# run experiments
	precision_res = {}
	relative_err_res = {}
	fnr_res = {}
	recall_res = {}
	
	for dataset_name, params in experiment.iteritems():
		logging.info('testing dataset ' + dataset_name)
		dataset = Dataset(datapath[dataset_name][0])
		dataset.mine_freq_item(datapath[dataset_name][1], float('inf'))
		
		for eps in eps_list:
			for k in params:
				#print eps, k, k1
				exp = ExperimentAlgo3(dataset)
				expresults = exp.run_experiment(eps, k, repeat=20)
				
				precision_res[(eps, dataset_name, k)] = [np.mean(expresults['precision']), stderr(expresults['precision'])]
				relative_err_res[(eps, dataset_name, k)] = [np.mean(expresults['relative_err']), stderr(expresults['relative_err'])]
				fnr_res[(eps, dataset_name, k)] = [np.mean(expresults['fnr']), stderr(expresults['fnr'])]
				recall_res[(eps, dataset_name, k)] = [np.mean(expresults['recall']), stderr(expresults['recall'])]

	# print results
	f = sys.stdout #open('method1_precision.txt', 'w')
	print '\n PRECISION \n\n'
	f.write('eps\t')
	for dataset_name, param in sorted(experiment.iteritems()):
		for k in param:
			f.write('(%s,%d)\tse\t' % (dataset_name, k))
	f.write('\n')
	
	for eps in eps_list:
		f.write('%.3f\t'%eps)
		for dataset_name, params in sorted(experiment.iteritems()):
			for k in params:
				f.write('%f\t%f\t' % (precision_res[(eps, dataset_name, k)][0], precision_res[(eps, dataset_name, k)][1]))
		f.write('\n')
	
	f = sys.stdout #open('method1_re.txt', 'w')
	print '\n RELATIVE ERR \n\n'
	f.write('eps\t')
	for dataset_name, param in sorted(experiment.iteritems()):
		for k in param:
			f.write('(%s,%d)\tse\t' % (dataset_name, k))
	f.write('\n')

	for eps in eps_list:
		f.write('%.3f\t'%eps)
		for dataset_name, params in sorted(experiment.iteritems()):
			for k in params:
				f.write('%f\t%f\t' % (relative_err_res[(eps, dataset_name, k)][0], relative_err_res[(eps, dataset_name, k)][1]))
		f.write('\n')
		
	f = sys.stdout #open('method1_fnr.txt', 'w')
	print '\n FNR \n\n'
	f.write('eps\t')
	for dataset_name, param in sorted(experiment.iteritems()):
		for k in param:
			f.write('(%s,%d)\tse\t' % (dataset_name, k))
	f.write('\n')
	
	for eps in eps_list:
		f.write('%.3f\t'%eps)
		for dataset_name, params in sorted(experiment.iteritems()):
			for k in params:
				f.write('%f\t%f\t' % (fnr_res[(eps, dataset_name, k)][0], fnr_res[(eps, dataset_name, k)][1]))
		f.write('\n')
	
	f = sys.stdout #open('method1_recall.txt', 'w')
	print '\n RECALL \n\n'
	f.write('eps\t')
	for dataset_name, param in sorted(experiment.iteritems()):
		for k in param:
			f.write('(%s,%d)\tse\t' % (dataset_name, k))
	f.write('\n')

	for eps in eps_list:
		f.write('%.3f\t'%eps)
		for dataset_name, params in sorted(experiment.iteritems()):
			for k in params:
				f.write('%f\t%f\t' % (recall_res[(eps, dataset_name, k)][0], recall_res[(eps, dataset_name, k)][1]))
		f.write('\n')
	
	
	
if __name__ == '__main__':
	main()


