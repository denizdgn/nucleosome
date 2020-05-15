
import numpy as np
import seaborn as sns
import pandas as pd
import os

import matplotlib.pyplot as plt
#pathx="/Volumes/DENIZ/kalelab/JMB/Yawen/cenpa_601_AB/analysisNuc/10_mutual_information"
#pathx="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_1/h3/analysisNuc/10_mutual_information"

pathx="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_1/cenpa/analysisNuc/10_mutual_information"
mini = np.load(f'{pathx}/mi.npy')

plt.close("all")
sns.heatmap(mini, cmap='inferno_r')
plt.xlabel('Base pair')
plt.ylabel('Base pair')
#plt.title('Normalized Mutual Information')
plt.title('Mutual Information')
plt.show()


path_write="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_1/cenpa/analysisNuc/00_graph/correlation"
plt.savefig(f'{path_write}/mi.png', bbox_inches='tight', dpi=300)
#plt.show()



#np.savetxt("/Users/denizdogan/Desktop/yy/10_mutual_information/foo.csv", mini, delimiter=",")


pathx="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_1/cenpa/analysisNuc/10_mutual_information"
mi = np.load(f'{pathx}/mi.npy')

plt.close("all")
sns.heatmap(mini, cmap='inferno_r')
plt.xlabel('Base pair')
plt.ylabel('Base pair')
#plt.title('Normalized Mutual Information')
plt.title('Mutual Information')
plt.show()



# -- Inferno -- #

pathx="/Volumes/DENIZ/kalelab/JMB/Yawen/cenpa_601_AB/analysisNuc/10_mutual_information"


for round in ["round_1","round_2","round_3"]:
    pathx=f'/Volumes/DENIZ/kalelab/Project_biowulf/MDx/{round}/h3/analysisNuc'
    mi = np.load(f'{pathx}/10_mutual_information/mi.npy')
    plt.close("all")
    sns.set_context('poster')
    ticklabels = ['SHL–7', 'SHL–6', 'SHL–5', 'SHL–4', 'SHL–3', 'SHL–2', 'SHL–1', 'Dyad', 'SHL+1', 'SHL+2', 'SHL+3', 'SHL+4', 'SHL+5', 'SHL+6', 'SHL+7']
    shl = np.linspace(0, 145, 15)
    #sns.heatmap(mi, cmap='tab20c', cbar_kws={'label': "Mutual Information"})
    sns.heatmap(mi, cmap='inferno_r', cbar_kws={'label': "Mutual Information"})
    plt.xticks(shl, ticklabels)
    plt.yticks(shl, ticklabels)
    #plt.show()
    os.makedirs(f'{pathx}/00_graph/mutual_info/', exist_ok=True)
    plt.savefig(f'{pathx}/00_graph/mutual_info/mi_2.png', bbox_inches='tight', dpi=300)





