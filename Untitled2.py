#!/usr/bin/env python
# coding: utf-8

# In[2]:


from Bio import SeqIO
from collections import Counter
import itertools
from itertools import permutations    #连续返回由 iterable 元素生成长度为 r 的排列。
import time
import pandas as pd
from math import sqrt
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import squiggle

# from Bio import SeqIO
# from collections import Counter
# import itertools
# from itertools import permutations
# from math import sqrt

start=time.perf_counter()
# 1) Calculate degree/score of twenty amino acids (R)
def sum_digits(your_string):
    return sum(int(x) for x in your_string if '0' <= x <= '9')


P = [''.join(P) for P in permutations('1234')]
# print(P)
N = ['A', 'T', 'C', 'G']
A = ['GATGAC', 'GAAGAG', 'CGTCGCCGACGGAGAAGG', 'CATCAC', 'AAAAAG', 'TATTAC', 'TTTTTC',
     'TGG', 'ATTATCATA', 'TTATTGCTTCTCCTACTG', 'GTTGTCGTAGTG', 'GCTGCCGCAGCG',

     'GGTGGCGGAGGG', 'CCTCCCCCACCG', 'ATG', 'TGTTGC', 'TCTTCCTCATCGAGTAGC', 'ACTACCACAACG', 'CAACAG',
     'AATAAC']
C = [2, 2, 6, 2, 2, 2, 2, 1, 3, 6, 4, 4, 4, 4, 1, 2, 6, 4, 2, 2]
R = [[0 for x in range(len(A))] for y in range(24)]
for loop in range(24):
    for i in range(20):
        st = A[i]
        for j in range(4):
            st = st.replace(N[j], P[loop][j])
        R[loop][i] = round(sum_digits(st) / C[i], 3)
    # print('', R[loop][:])

# 2) Generate feature vector (FV)
fnam = str("all_virus2.fasta")
sn = 2
AA = ['D','E','R','H','K','Y','F','W','I','L','V','A','G','P','M','C','S','T','Q','N']

FV = []
TN = []
i = -1
for record in SeqIO.parse(r'F:\testprogram\venv\all_virus2.fasta', "fasta"):
    i = i + 1
    fv = []
    TN.append(record.id)
    #统计频率
    path = r'D:/3333.xlsx'
    df = pd.read_excel(path, index_col=0)
    dfList = df.values
    # print(dfList)
    data = pd.DataFrame(df)
    # print(data)
    # 画频率图

# plt.axis()
# plt.ylabel('Frequencey')
# xlable=['D', 'E', 'R', 'H', 'K', 'Y', 'F', 'W', 'I', 'L', 'V', 'A', 'G', 'P', 'M', 'C', 'S', 'T', 'Q', 'N']
# # plt.rcParams['font.sans-serif'] = ['SimHei']
#
# plt.xticks(rotation=30)
# plt.plot(xlable,data)
# plt.legend([1,2,3,4,5,6,7])
# plt.show()


    # for i in range(7):
    for j in range(20):
        # print(len(record.seq)
        fv.append((round((df.iloc[i][j]* R[sn - 1][j])*100 / (len(record.seq)-2),4)))
        # print(df.iloc[i][j]* R[sn - 1][j]/ ((len(record.seq)-2)))

    FV.append(fv)
    print(FV)
n = i + 1


# 3) Calculate pairwise distance matrix
f2 = open('M-%s.tsv' % fnam, 'w')
f2.write('%s' % n)
SDM = [[0 for x in range(n)] for y in range(n)]
for i in range(n):
    f2.write('\n%-12s' % TN[i])
    for j in range(n):
        s = 0
        for k in range(20):
            s = s + pow(FV[i][k] - FV[j][k], 2)
        f2.write('\t%8.4f' % sqrt(s))
f2.close()
end=time.perf_counter()

print(end-start)


# In[ ]:




