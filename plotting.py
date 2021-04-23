import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
# plt.figure(figsize=(8, 4))
# fig, ax = plt.subplots(1,2)
# sns.heatmap(corr_omega, ax=ax[0])
# sns.heatmap(corr_domega, ax=ax[1]);

f = plt.figure(figsize=(10,10))
ax = f.add_subplot(221)
ax1 = f.add_subplot(222)
ax2 = f.add_subplot(223)
ax3 = f.add_subplot(224)
sns.heatmap(cov_matrix, ax=ax)
sns.heatmap(corr_omega, ax=ax1)
sns.heatmap(corr_domega, ax=ax2)
sns.heatmap(corr_any, ax=ax3)

ax1.title.set_text('Cross-node AFV')
ax2.title.set_text('Cross-node RoCoF')
ax3.title.set_text('Cross-node AV')
ax.title.set_text('Covariance matrix of disturbances')
plt.savefig('corr_matrices.png')




plt.plot(can_alpha, alpha_rates[:,0], 'r',label = 'AFV')
plt.plot(can_alpha, alpha_rates[:,1], 'b', label = 'RoCoF')
plt.plot(can_alpha, alpha_rates[:,2], 'g', label= 'AV')

