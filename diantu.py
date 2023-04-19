#!/usr/bin/env python
# coding: utf-8

# In[84]:


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.pyplot import figure

figure(figsize=(5,5), dpi=120)

from scipy.stats import t

# plt.close('all')


data= np.array([[19.5646,20.8689,37.5014,17.3456,22.6196,26.76,43.2661,20.3806,31.7648,98.8654,63.2253,33.2029,38.3817,19.3037,19.3973,44.639,55.7088,37.2111,20.9876,20.7083],[20.9503,21.4041,41.7835,19.1351,20.0124,25.3235,39.4299,20.5957,29.4464,96.2636,59.7482,40.9442,42.4603,21.3032,20.8679,46.3679,59.2709,36.5878,21.6495,17.9855],[19.499,17.1581,41.0018,16.8925,16.7347,27.3666,41.2342,21.3667,31.2182,102.6298,66.494,42.0527,43.0222,25.1436,19.0457,47.5047,63.4525,32.2509,18.3717,16.0474],[24.3314,27.8525,51.5257,18.8875,24.4263,20.6757,36.0948,15.1448,40.1646,77.8879,41.2935,34.5703,46.4156,30.0364,13.8419,24.9064,66.8618,35.686,26.634,22.7409],[19.9477,20.8522,37.4523,17.1155,22.2942,27.3473,42.7093,20.5607,31.7355,99.7853,62.6269,33.0,38.5489,18.9864,19.5357,44.0425,56.3009,37.5289,20.9108,20.5758],[20.0084,21.7701,40.2047,17.8473,21.3372,26.1309,42.5856,20.3389,31.2685,97.8862,63.1527,34.9698,40.6258,19.5218,19.2483,45.2366,56.9547,36.7466,20.6527,19.6309],[20.2825,21.6196,40.1894,17.9263,21.2361,26.2143,42.257,20.0538,31.5843,98.5646,62.7901,34.8234,39.9798,19.5022,18.9976,44.7999,57.0888,36.9627,20.4053,19.6165]]
) 
data[:,0] = (data[:,0] - np.min(data[:,0]))/(np.max(data[:,0]) - np.min(data[:,0]))
data[:,1] = (data[:,1] - np.min(data[:,1]))/(np.max(data[:,1]) - np.min(data[:,1]))

class_colours = ["r",  "r", "g", "g", "b", "b", "b"]
custom_annotations=['NC_045512.2','NC_004718.3','NC_019843.3','NC_002549.1','MN996532.2 ','MG772933.1','MG772934.1']


for i, point in enumerate(data): 
    plt.scatter(point[0], point[1],label=custom_annotations[i],edgecolors='black', linewidths=1, alpha=0.75)
    
    #     plt.annotate(i+1, (data[i,0], data[i,1]))



# plt.axvline(0, c=(.5, .5, .5), ls= '--')
# plt.axhline(0, c=(.5, .5, .5), ls= '--')

plt.legend(ncol=1)


plt.savefig('D:dian.png',dpi=600)
plt.show()


# In[3]:


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.pyplot import figure


plt.close('all')


data= np.array([[19.5646,20.8689,37.5014,17.3456,22.6196,26.76,43.2661,20.3806,31.7648,98.8654,63.2253,33.2029,38.3817,19.3037,19.3973,44.639,55.7088,37.2111,20.9876,20.7083],[20.9503,21.4041,41.7835,19.1351,20.0124,25.3235,39.4299,20.5957,29.4464,96.2636,59.7482,40.9442,42.4603,21.3032,20.8679,46.3679,59.2709,36.5878,21.6495,17.9855],[19.499,17.1581,41.0018,16.8925,16.7347,27.3666,41.2342,21.3667,31.2182,102.6298,66.494,42.0527,43.0222,25.1436,19.0457,47.5047,63.4525,32.2509,18.3717,16.0474],[24.3314,27.8525,51.5257,18.8875,24.4263,20.6757,36.0948,15.1448,40.1646,77.8879,41.2935,34.5703,46.4156,30.0364,13.8419,24.9064,66.8618,35.686,26.634,22.7409],[19.9477,20.8522,37.4523,17.1155,22.2942,27.3473,42.7093,20.5607,31.7355,99.7853,62.6269,33.0,38.5489,18.9864,19.5357,44.0425,56.3009,37.5289,20.9108,20.5758],[20.0084,21.7701,40.2047,17.8473,21.3372,26.1309,42.5856,20.3389,31.2685,97.8862,63.1527,34.9698,40.6258,19.5218,19.2483,45.2366,56.9547,36.7466,20.6527,19.6309],[20.2825,21.6196,40.1894,17.9263,21.2361,26.2143,42.257,20.0538,31.5843,98.5646,62.7901,34.8234,39.9798,19.5022,18.9976,44.7999,57.0888,36.9627,20.4053,19.6165]]) 

custom_annotations = ["K464E", "K472E", "R470E", "K464A", "M155E", "K472A", "M155A"]
class_colours = ["a", "b", "c", "d", "e", "f", "g"]
class_desc = ['description 1', 'description 2', 'description 3']
dataclass = np.array([0,0,0,0,0,0,1])
colorclass = ['r','g','b']

for datac in range(3):
    plt.scatter(data[dataclass==datac][:,0],data[dataclass==datac][:,1], c=colorclass[datac], label=class_desc[datac])

for i, point in enumerate(data): 
#     plt.scatter(point[0], point[1], marker='o', label=custom_annotations[i], c=class_colours[i], edgecolors='black', linewidths=1, alpha=0.75)
     plt.annotate(custom_annotations[i], (data[i,0], data[i,1]))



plt.axvline(0, c=(.5, .5, .5), ls= '--')
plt.axhline(0, c=(.5, .5, .5), ls= '--')

plt.legend(fontsize=12)

plt.show()


# In[ ]:




