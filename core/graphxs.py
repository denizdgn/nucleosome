#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.legend_handler import HandlerPatch
from collections import OrderedDict
from matplotlib.colors import ListedColormap
import seaborn as sns
import pandas as pd


class contact():
    def heatmaps(filex,selection1,selection2):

    #path="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_2/h3/analysisNuc/02_contact"
    #pathr="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/"
    #os.chdir(path)

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        path=os.getcwd()
        if not os.path.exists("00_graph/contact/heatmap"):
            os.makedirs("00_graph/contact/heatmap", exist_ok=True)
        pathx = path + "/00_graph/contact/heatmap"


        df= pd.read_csv(filex)

        ind = []
        for j in range(len(df)):
            if all(df.iloc[j, 1:] == 0):
                ind.append(j)

        df.drop(df.index[ind], inplace=True)

        df=df.reset_index(drop=True)

        l=list(range(1,len(df.columns)))

        df_ultimate=pd.DataFrame(columns=["resid", "number", "time"])
        for i in range(len(df)):
            lf_2 = [[df["resid"][i].tolist()] * len(l), df.iloc[i, 1:-1].tolist(), l]
            df_2=pd.DataFrame(lf_2).T
            df_2.columns = ["resid", "number", "time"]
            df_ultimate = pd.concat([df_ultimate,df_2])

        #print(j)


        df_ultimate.resid = df_ultimate.resid.astype(int)
        df_ultimate.time = df_ultimate.time.astype(int)


        heatmap1_data = pd.pivot_table(df_ultimate, values='number',
                             index=['resid'],
                             columns='time')

        cmap = sns.cubehelix_palette(n_colors=100,dark=0,light=0.975,gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        #flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
        #cmap = ListedColormap(sns.color_palette(flatui).as_hex())

        #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [ "#f9f5f7","#143943"])

        fig,ax=plt.subplots(1, 1,figsize=(25,25))
        #sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
        #            ax=ax)

        sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                    ax=ax, vmin=0, vmax=150)

        for _, spine in ax.spines.items():
            spine.set_visible(True)

        for t in ax.xaxis.get_ticklines(): t.set_color('silver')
        for t in ax.yaxis.get_ticklines(): t.set_color('silver')

        plt.title(f'Contact Number Per Residue - {selection1} to {selection2}',
                  fontdict=dict(weight='bold'),fontsize = 24)
        plt.ylabel("Residue Number", fontdict=dict(weight='bold'),fontsize = 20)
        plt.xlabel("Frame Number", fontdict=dict(weight='bold'),fontsize = 20)
        yticks=ax.get_yticklabels()
        ax.set_yticklabels(yticks,fontdict=dict(weight='bold',fontsize=14),rotation=0)

        xticks_max=max(df_ultimate.time)
        #xticks_min=min(df_ultimate.time)
        xticks=np.arange(0, xticks_max, 100)
        ax.xaxis.set_ticks(xticks)
        ax.set_xticklabels(xticks,fontdict=dict(fontsize=14))
        file_name = f'{sel1}_to_{sel2}'
        plt.savefig(f'{pathx}/{file_name}.png',bbox_inches='tight',dpi=300)

        plt.close()
