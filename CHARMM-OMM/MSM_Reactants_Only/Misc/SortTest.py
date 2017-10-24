from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

phenos = np.array([128, 20, 0, 144, 4, 16, 160, 136, 192, 52, 128, 20, 0, 4, 16, 144, 130, 136, 132, 22,
128, 160, 4, 0, 32, 36, 132, 136, 164, 130, 128, 22, 4, 0, 144, 160, 54, 130, 178, 132,
128, 4, 0, 136, 132, 68, 196, 130, 192, 8, 128, 4, 0, 20, 22, 132, 144, 192, 130, 2,
128, 4, 0, 132, 20, 136, 144, 192, 64, 130, 128, 4, 0, 144, 132, 28, 192, 20, 16, 136,
128, 6, 4, 134, 0, 130, 160, 132, 192, 2,  128, 4, 0, 132, 68, 160, 192, 36, 64,
128, 4, 0, 136, 192, 8, 160, 12, 36, 128, 4, 0, 22, 20, 144, 86, 132, 82, 160,
128, 4, 0, 132, 20, 192, 144, 160, 68, 64, 128, 4, 0, 132, 160, 144, 136, 192, 68, 20])

values = np.sort(phenos)
indexes = np.arange(len(values))
width = 1
plt.bar(indexes, values, width)
# plt.xticks(indexes + width * 0.5, labels)
plt.savefig('Test')