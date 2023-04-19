
from Bio import SeqIO
from collections import Counter
import itertools
from itertools import permutations    #连续返回由 iterable 元素生成长度为 r 的排列。
from math import sqrt
import time

covid = SeqIO.read(r'D:/data\data3\COVID-19.fasta', 'fasta')
sars = SeqIO.read(r'D:/data\data3\SARS.fasta', 'fasta')
mers = SeqIO.read(r'D:/data\data3\MERS.fasta', 'fasta')
ebola = SeqIO.read(r'D:/data\data3\Ebola.fasta', 'fasta')
RaTG13= SeqIO.read(r'D:data\data3\RaTG13.fasta', 'fasta')
CoVZC45= SeqIO.read(r'D:\data\data3\CoVZC45.fasta', 'fasta')
CoVZXC21= SeqIO.read(r'D:\data\data3\CoVZXC21.fasta', 'fasta')
## 合并文件
all_virus = [covid, sars, mers, ebola,RaTG13,CoVZC45,CoVZXC21]
def getKmers(sequence, size=3):
    return [sequence[x:x+size] for x in range(len(sequence) - size + 1)]
# # function to convert sequence strings into k-mer words, default size = 6 (hexamer words)
fa_seq=SeqIO.parse(r'D:/data/data3/all_virus2.fasta', 'fasta')
t=0

for fa in SeqIO.parse("D:/data/data3/all_virus2.fasta", "fasta"):
    seqs=[str(fa.seq)]
    # print(seqs)
   # print(len(seqs[0]))
    kmermartix = getKmers(seqs[0])
   # print(kmermartix)
    result=Counter(kmermartix)#利用collection包下得Counter的类
    # print(result)








from Bio import SeqIO
from collections import Counter
import itertools
from itertools import permutations    #连续返回由 iterable 元素生成长度为 r 的排列。
from math import sqrt
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import squiggle
# 1) Calculate degree/score of twenty amino acids (R)
def sum_digits(your_string):
    return sum(int(x) for x in your_string if '0' <= x <= '9')
start=time.time()

P = [''.join(P) for P in permutations('1234')]

N = ['A', 'T', 'C', 'G']
A = list(result.keys())
resultValue=list(result.values())
print('values',resultValue)

C = [2, 2, 6, 2, 2, 2, 2, 1, 3, 6, 4, 4, 4, 4, 1, 2, 6, 4, 2, 2]
from Bio import SeqIO
from collections import Counter
import itertools
from itertools import permutations    #连续返回由 iterable 元素生成长度为 r 的排列。
from math import sqrt
import time
#sample1
# covid = SeqIO.read(r'D:/data\data3\COVID-19.fasta', 'fasta')
# sars = SeqIO.read(r'D:/data\data3\SARS.fasta', 'fasta')
# mers = SeqIO.read(r'D:/data\data3\MERS.fasta', 'fasta')
# ebola = SeqIO.read(r'D:/data\data3\Ebola.fasta', 'fasta')
# RaTG13= SeqIO.read(r'D:data\data3\RaTG13.fasta', 'fasta')
# CoVZC45= SeqIO.read(r'D:\data\data3\CoVZC45.fasta', 'fasta')
# CoVZXC21= SeqIO.read(r'D:\data\data3\CoVZXC21.fasta', 'fasta')
# ## 合并文件
# all_virus = [covid, sars, mers, ebola,RaTG13,CoVZC45,CoVZXC21]


def getKmers(sequence, size=3):
    return [sequence[x:x+size] for x in range(len(sequence) - size + 1)]
# # function to convert sequence strings into k-mer words, default size = 6 (hexamer words)
fa_seq=SeqIO.parse(r'D:/data/data3/all_virus2.fasta', 'fasta')
t=0
from collections import Counter
for fa in SeqIO.parse("D:/data/data3/all_virus2.fasta", "fasta"):
    seqs=[str(fa.seq)]
    # print(seqs)
   # print(len(seqs[0]))
    kmermartix = getKmers(seqs[0])
   # print(kmermartix)
    result=Counter(kmermartix)#利用collection包下得Counter的类
    # print(result)  #创建一个[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
R=[[0 for x in range(len(A))] for y in range(24)]
#print(R)
for loop in range(24):
    for i in range(20):#20个氨基酸？
        st = A[i]
        for j in range(4):
            st = st.replace(N[j], P[loop][j])#atcg权重赋值
        R[loop][i] = round(sum_digits(st) / C[i], 3)
    # print('', R[loop][:21])



# 2) Generate feature vector (FV)

fnam = str("all_virus2.fasta")
# sn = int(input('Enter Scheme No.: '))
sn=1#权值分配的方案
AA = ['D', 'E', 'R', 'H', 'K', 'Y', 'F', 'W', 'I', 'L', 'V', 'A', 'G', 'P', 'M', 'C', 'S', 'T', 'Q', 'N']
FV = []
TN = []
i = -1
#
# def nt_stat(s):
#     return s.count("A"), s.count("G"), s.count("C"), s.count("T")
for record in SeqIO.parse(r'F:\testprogram\venv\all_virus2.fasta', "fasta"):
    i = i + 1
    # print(nt_stat(record))
    fv = []
    TN.append(record.id)
    for j in range(20):
        # fv.append(round((((record.seq).count(AA[j])) * R[sn - 1][j] * 100) / (len(record.seq)), 4))
        print((record.seq).count(AA[j]))

    FV.append(fv)
    print(FV)
n = i + 1

# print(n)


#3) Calculate pairwise distance matrix
f2=open('DM-%s.tsv'%fnam,'w')
f2.write('%s'%n)
# SDM=[[0 for x in range(n)] for y in range(n)]
# import numpy as np
#
#
# def dtw_distance(ts_a, ts_b, d=lambda x, y: abs(x - y), mww=10000):
#     """Computes dtw distance between two time series
#
#     Args:
#         ts_a: time series a
#         ts_b: time series b
#         d: distance function
#         mww: max warping window, int, optional (default = infinity)
#
#     Returns:
#         dtw distance
#     """
#
#     # Create cost matrix via broadcasting with large int
#     ts_a, ts_b = np.array(ts_a), np.array(ts_b)
#     M, N = len(ts_a), len(ts_b)
#     cost = np.ones((M, N))
#
#     # Initialize the first row and column
#     cost[0, 0] = d(ts_a[0], ts_b[0])
#     for i in range(1, M):
#         cost[i, 0] = cost[i - 1, 0] + d(ts_a[i], ts_b[0])
#
#     for j in range(1, N):
#         cost[0, j] = cost[0, j - 1] + d(ts_a[0], ts_b[j])
#
#     # Populate rest of cost matrix within window
#     for i in range(1, M):
#         for j in range(max(1, i - mww), min(N, i + mww)):
#             choices = cost[i - 1, j - 1], cost[i, j - 1], cost[i - 1, j]
#             cost[i, j] = min(choices) + d(ts_a[i], ts_b[j])
#
#     # Return DTW distance given window
#     return cost[-1, -1]
# for i in range(n):
#     f2.write('\n%-12s'%TN[i])
#     for j in range(n):
#         s=0
#         for k in range(20):
#             s=s+dtw_distance(FV[i][k]，FV[j][k]）
#         f2.write('\t%8.4f'%sqrt(s))
# # f2.close()
# end=time.time()
# print(end-start)
