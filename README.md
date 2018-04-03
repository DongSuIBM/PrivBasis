# PrivBasis: Frequent Itemset Mining with Differential Privacy
PrivBasis is a differentially private frequent itemset mining algorithm, which leverages a novel notion called basis sets. A theta-basisset has the property that any itemset with frequency higher than theta is a subset of some basis.  PrivBasis first privately constructs a basis set and then uses it to find the most frequent itemsets.  

For more details, please see our paper:
Ninghui Li, Wahbeh H. Qardaji, Dong Su, Jianneng Cao: [PrivBasis: Frequent Itemset Mining with Differential Privacy](http://vldb.org/pvldb/vol5/p1340_ninghuili_vldb2012.pdf). PVLDB 5(11): 1340-1351 (2012).  

## Assumptions
- The dataset is assumed to be space separated.
- The max clique size (or basis length) is less than 12 as we claimed in the Section 4.2. 

## How to run
The entrance of the experiment is the PrivBasisMain.py file. Please go to line 300 to change the dataset path, minimum support count, repeat time and k values. When the experiment is done, the evaluation result will be printed in the terminal. 

## Environment requirements
- Python 2.7.12
- Numpy 1.12.1
- Pandas 0.19.1

Maintainer:
Dong Su, <sudong.tom@gmail.com>