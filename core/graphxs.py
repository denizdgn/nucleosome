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
import matplotlib.cm as cm
import re


class nuc_rc_sm():
    def heatmaps_nuc_e(filex):

        path = os.getcwd()
        if not os.path.exists("00_graph/nuc_e/heatmap"):
            os.makedirs("00_graph/nuc_e/heatmap", exist_ok=True)
        pathx = path + "/00_graph/nuc_e/heatmap"

        df = pd.read_csv(filex)

        abs_path = os.path.abspath(filex)
        rox = abs_path.split("MDx/")[1].split(f'/analysisNuc')[0]
        ro = rox.split("/")[0]
        cx = rox.split("/")[1]

        heatmap1_data = df.T
        heatmap1_data.drop('time', axis=0, inplace=True)

        cmap = sns.cubehelix_palette(n_colors=1800, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        # flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
        # cmap = ListedColormap(sns.color_palette(flatui).as_hex())

        # cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [ "#f9f5f7","#143943"])

        fig, ax = plt.subplots(1, 1, figsize=(25, 25))
        # sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
        #            ax=ax)

        sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                    ax=ax, vmin=0, vmax=1800)

        for _, spine in ax.spines.items():
            spine.set_visible(True)

        for t in ax.xaxis.get_ticklines(): t.set_color('silver')
        for t in ax.yaxis.get_ticklines(): t.set_color('silver')

        plt.title(f'Bending Energy ({cx})',
                  fontdict=dict(weight='bold'), fontsize=24)
        plt.ylabel("", fontdict=dict(weight='bold'), fontsize=20)
        plt.xlabel("Time Scale (ns)", fontdict=dict(weight='bold'), fontsize=20)
        yticks = ax.get_yticklabels()
        ax.set_yticklabels(yticks, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        xticks = np.arange(0, 5002, 500)
        ax.xaxis.set_ticks(xticks)
        xtick_label = np.arange(0, 502, 50)
        ax.set_xticklabels(xtick_label, fontdict=dict(fontsize=14))
        file_name = f'{cx}_{ro}_nuc_e'
        # plt.show()

        plt.savefig(f'{pathx}/{file_name}.png', bbox_inches='tight', dpi=300)

        plt.close()

    def heatmaps_rc_sm(filex):

        path = os.getcwd()
        if not os.path.exists("00_graph/nuc_rc_sm/heatmap"):
            os.makedirs("00_graph/nuc_rc_sm/heatmap", exist_ok=True)
        pathx = path + "/00_graph/nuc_rc_sm/heatmap"

        df = pd.read_csv(filex)
        df.drop('mean', axis=1, inplace=True)

        abs_path = os.path.abspath(filex)
        rox = abs_path.split("MDx/")[1].split(f'/analysisNuc')[0]
        ro = rox.split("/")[0]
        cx = rox.split("/")[1]
        sm = abs_path.split("08_nuc_rc_sm_merged/")[1].split("sm")[1].split(".cs")[0]

        heatmap1_data = df.T
        heatmap1_data.drop('time', axis=0, inplace=True)

        cmap = sns.cubehelix_palette(n_colors=80, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        fig, ax = plt.subplots(1, 1, figsize=(35, 35))

        sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                    ax=ax, vmin=0, vmax=80)

        for _, spine in ax.spines.items():
            spine.set_visible(True)

        for t in ax.xaxis.get_ticklines(): t.set_color('silver')
        for t in ax.yaxis.get_ticklines(): t.set_color('silver')

        plt.title(f'Radius of Curvature ({cx} - sm{sm})',
                  fontdict=dict(weight='bold'), fontsize=26)
        plt.ylabel("", fontdict=dict(weight='bold'), fontsize=20)
        plt.xlabel("Time Scale (ns)", fontdict=dict(weight='bold'), fontsize=24)
        yticks = ax.get_yticklabels()
        ax.set_yticklabels(yticks, fontdict=dict(weight='bold', fontsize=16), rotation=0)

        xticks = np.arange(0, 5002, 500)
        ax.xaxis.set_ticks(xticks)
        xtick_label = np.arange(0, 502, 50)
        ax.set_xticklabels(xtick_label, fontdict=dict(fontsize=14))
        file_name = f'{cx}_{ro}_sm{sm}'
        # plt.show()

        plt.savefig(f'{pathx}/{file_name}.png', bbox_inches='tight', dpi=300)

        plt.close()

class contact():
    def multiple_heatmaps(pathq, selection1, selection2):
        pathv=os.getcwd()
        pathb = pathv+"/"+pathq

        print(pathb)
        roundsl = {"round_1": "601", "round_2": "Alpha Satellite Sequence", "round_3": "601L"}
        cenpa_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L",
                        "chainE": "CENPA.R",
                        "chainF": "H4.R", "chainG": "H2A.R",
                        "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                        "nucleic": "Nucleic"}
        cenpa_h2az_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.Z.L", "chainD": "H2B.L",
                             "chainE": "CENPA.R",
                             "chainF": "H4.R", "chainG": "H2A.Z.R",
                             "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                             "nucleic": "Nucleic"}
        h3_chains = {"chainA": "H3.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L", "chainE": "H3.R",
                     "chainF": "H4.R", "chainG": "H2A.R",
                     "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                     "nucleic": "Nucleic"}

        sel1x = selection1.replace(" ", "")
        sel2x = selection2.replace(" ", "")

        file_name = f'{sel1x}_to_{sel2x}'

        if not os.path.exists("00_graph/contact/heatmap"):
            os.makedirs("00_graph/contact/heatmap", exist_ok=True)
        pathx = pathb + "/00_graph/contact/heatmap"

        round = pathb.split("MDx/")[1]
        cx = ["h3", "cenpa", "cenpa-h2az"]
        cx_rmsf = OrderedDict()
        sells = OrderedDict()

        for cxy in cx:

            filex = f'/Volumes/DENIZ/kalelab/Project_biowulf/MDx/{round}/{cxy}/analysisNuc/06_contact_merged/{file_name}/{file_name}.csv'

            if cxy == "cenpa":
                sel1 = cenpa_chains[sel1x]
                sel2 = cenpa_chains[sel2x]
                print(sel1, sel2)
            elif cxy == "cenpa-h2az":
                sel1 = cenpa_h2az_chains[sel1x]
                sel2 = cenpa_h2az_chains[sel2x]
                print(sel1, sel2)
            elif cxy == "h3":
                sel1 = h3_chains[sel1x]
                sel2 = h3_chains[sel2x]
                print(sel1, sel2)

            sells[cxy] = [sel1,sel2]

            df = pd.read_csv(filex)
            pathh = os.path.abspath(filex)
            ro = pathh.split("MDx/")[1].split(f'/{cxy}')[0]
            file = os.path.basename(f'{filex}')
            chain = file.split(".csv")[0].split("chain")[1].split("_")[0]

            time = []
            for i in range(1, len(df.columns) - 1):
                time.append(df.columns[i].split("number_")[1])
            print(cxy, ro, chain)

            if chain == "J":
                df = df[::-1]
                df.reset_index(drop=True, inplace=True)
            if chain == "J" and cxy == "h3" and ro == "round_2":
                df.resid = df.resid - 73
            if chain == "I" and cxy == "h3" and ro == "round_2":
                df.resid = df.resid - 73
            if chain == "I" and cxy == "cenpa-h2az" and ro == "round_1":
                df.resid = df.resid - 73
            if chain == "J" and cxy == "cenpa-h2az" and ro == "round_1":
                df.resid = df.resid - 219
            if chain == "I" and cxy == "cenpa" and ro == "round_2":
                df.resid = df.resid - 73
            if chain == "J" and cxy == "cenpa" and ro == "round_2":
                df.resid = df.resid - 73
            if chain == "I" and cxy == "cenpa-h2az" and ro == "round_2":
                df.resid = df.resid - 73
            if chain == "J" and cxy == "cenpa-h2az" and ro == "round_2":
                df.resid = df.resid - 73
            if chain == "C":
                if cxy == "cenpa-h2az":
                    df.resid = df.resid - 800

            df_ultimate = pd.DataFrame(columns=["resid", "number", "time"])
            for i in range(len(df)):
                lf_2 = [[df["resid"][i]] * len(time), df.iloc[i, 1:-1].tolist(), time]
                df_2 = pd.DataFrame(lf_2).T
                df_2.columns = ["resid", "number", "time"]
                df_2.time = [int(i) / 10 for i in df_2.time]
                df_ultimate = pd.concat([df_ultimate, df_2])
            df_ultimate.resid = df_ultimate.resid.astype(int)
            df_ultimate.time = df_ultimate.time.astype(int)
            df_ultimate.number = df_ultimate.number.astype(int)
            heatmap1_data = pd.pivot_table(df_ultimate, values='number',
                                           index=['resid'],
                                           columns='time')

            cx_rmsf[cxy] = heatmap1_data

        cmap = sns.cubehelix_palette(n_colors=80, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        fig, ax = plt.subplots(1, 3, figsize=(60, 25))

        for i, cxy in enumerate(cx):

            if i < 2:
                sns.heatmap(cx_rmsf[cxy], cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                            ax=ax[i], vmin=0, vmax=80, cbar=False)
            else:
                sns.heatmap(cx_rmsf[cxy], cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                            ax=ax[i], vmin=0, vmax=80)

            for _, spine in ax[i].spines.items():
                spine.set_visible(True)

            for t in ax[i].xaxis.get_ticklines(): t.set_color('silver')
            for t in ax[i].yaxis.get_ticklines(): t.set_color('silver')

            ax[i].set_title(f'{cxy.upper()} Complex - ({sells[cxy][0]} to {sells[cxy][1]})', fontdict=dict(weight='bold'), fontsize=26)

            ax[i].set_ylabel(f'Residue Number ({sells[cxy][0]})', fontdict=dict(weight='bold'), fontsize=20)
            ax[i].set_xlabel("Time Scale (ns)", fontdict=dict(weight='bold'), fontsize=20)

            yticks = ax[i].get_yticklabels()
            ax[i].set_yticklabels(yticks, fontdict=dict(weight='bold', fontsize=16), rotation=0)

            xticks = np.arange(0, 502, 50)
            ax[i].xaxis.set_ticks(xticks)
            ax[i].set_xticklabels(xticks, fontdict=dict(fontsize=14))
            plt.suptitle(f' Contact number at {roundsl[round]} Sequence', y=1.05, fontweight="bold", fontsize=30,
                         fontdict=dict(weight='bold'))

        # plt.show()
        fig.tight_layout()
        logging.info(f'writing to {pathx} ...')
        plt.savefig(f'{pathx}/{file_name}.png', bbox_inches='tight', dpi=300)
        plt.close()

    def heatmap_barplot(filex, selection1, selection2):
        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        path = os.getcwd()
        if not os.path.exists("00_graph/contact/heatmap_barplot"):
            os.makedirs("00_graph/contact/heatmap_barplot", exist_ok=True)
        pathx = path + "/00_graph/contact/heatmap_barplot"

        df = pd.read_csv(filex)
        abs_path=os.path.abspath(filex)
        rox = abs_path.split("MDx/")[1].split(f'/analysisNuc')[0]
        #rox = abs_path.split("JMB/")[1].split(f'/analysisNuc')[0]
        ro = rox.split("/")[0]
        cx = rox.split("/")[1]

        ind = []
        for j in range(len(df)):
            if all(df.iloc[j, 1:] == 0):
                ind.append(j)

        df.drop(df.index[ind], inplace=True)

        df = df.reset_index(drop=True)

        df.drop('mean', axis=1, inplace=True)

        mean_row = df.iloc[:, 1:-1].mean(axis=1)
        df["mean"] = mean_row
        #sorted(df["mean"])




        df.resid = df.resid.astype(int)
        time = []
        for i in range(1, len(df.columns) - 1):
            time.append(df.columns[i].split("number_")[1])

        df_ultimate = pd.DataFrame(columns=["resid", "number", "time"])
        for i in range(len(df)):
            lf_2 = [[df["resid"][i]] * len(time), df.iloc[i, 1:-1].tolist(), time]
            df_2 = pd.DataFrame(lf_2).T
            df_2.columns = ["resid", "number", "time"]
            df_2.time = [int(i) / 10 for i in df_2.time]
            df_ultimate = pd.concat([df_ultimate, df_2])

        df_ultimate.resid = df_ultimate.resid.astype(int)
        df_ultimate.time = df_ultimate.time.astype(int)
        df_ultimate.number = df_ultimate.number.astype(int)

        heatmap1_data = pd.pivot_table(df_ultimate, values='number',
                                       index=['resid'],
                                       columns='time')

        cmap = sns.cubehelix_palette(n_colors=150, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        #data = df.groupby("mean").size()
        #norm = plt.Normalize(0, 150)
        #pal = plt.cm.PuBuGn(df["mean"].apply(norm))
        #pal = sns.color_palette("PuBuGn_d", 6000)
        #rank = data.argsort().argsort()


        #cmap=sns.color_palette("PuBuGn", 150)
        fig, ax_box = plt.subplots(1, 3, figsize=(25, 25), gridspec_kw={"width_ratios": [2, 1, 0.025]})

        sns.heatmap(heatmap1_data, cmap=cmap,  yticklabels=True,
                    cbar=False, ax=ax_box[0], vmin=0, vmax=150)

        #np.array(pal[::-1])[rank]
        sns.barplot(y="resid", x="mean", data=df, orient='h',hue="mean",dodge=False,
                    palette=cmap.colors, ax=ax_box[1])
        ax_box[1].legend().remove()

        fig.colorbar(ax_box[0].get_children()[0], cax=ax_box[2], orientation="vertical")

        ax_box[0].set_ylabel('Residue Number', fontdict=dict(weight='bold'), fontsize=18)
        ax_box[0].set_xlabel('Time Scale (ns)', fontdict=dict(weight='bold'), fontsize=16)

        for _, spine in ax_box[0].spines.items():
            spine.set_visible(True)

        for t in ax_box[0].xaxis.get_ticklines(): t.set_color('silver')
        for t in ax_box[0].yaxis.get_ticklines(): t.set_color('silver')
        xticks = np.arange(0, 502, 50)
        ax_box[0].xaxis.set_ticks(xticks)
        ax_box[0].set_xticklabels(xticks, fontdict=dict(fontsize=14))

        xticks0 = ax_box[0].get_xticklabels()
        ax_box[0].set_xticklabels(xticks0, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        yticks0 = ax_box[0].get_yticklabels()
        ax_box[0].set_yticklabels(yticks0, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        yticks1 = ax_box[1].get_yticklabels()
        ax_box[1].set_yticklabels(yticks1, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        ax_box[1].set_xlabel('Average Contact Number', fontdict=dict(weight='bold'), fontsize=16)
        ax_box[1].set_ylabel('')

        xticks1 = np.arange(0, 150, 20)
        ax_box[1].set_xticks(xticks1)
        ax_box[1].set_xticklabels(xticks1, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        for t in ax_box[1].xaxis.get_ticklines(): t.set_color('silver')
        for t in ax_box[1].yaxis.get_ticklines(): t.set_color('silver')
        fig.suptitle(f'{cx}: {selection1} to {selection2}', fontdict=dict(weight='bold'), fontweight="bold", fontsize=24,
                     y=0.92)
        file_name = f'{sel1}_to_{sel2}'
        plt.savefig(f'{pathx}/{file_name}.png', bbox_inches='tight', dpi=300)
        plt.close()

    def heatmaps(filex,selection1,selection2):

    #path="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_2/h3/analysisNuc/02_contact"
    #pathr="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/"
    #os.chdir(path)

        sel1x = selection1.replace(" ", "")
        sel2x = selection2.replace(" ", "")

        cenpa_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L", "chainE": "CENPA.R",
                    "chainF": "H4.R", "chainG": "H2A.R",
                    "chainH": "H2B.R", "chainM":"ABM", "chainN":"ABN", "chainX":"ABX", "chainM":"ABY",
                        "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein", "nucleic": "Nucleic"}

        cenpa_h2az_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.Z.L", "chainD": "H2B.L", "chainE": "CENPA.R",
                    "chainF": "H4.R", "chainG": "H2A.Z.R",
                    "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein", "nucleic": "Nucleic"}

        h3_chains = {"chainA": "H3.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L", "chainE": "H3.R",
                    "chainF": "H4.R", "chainG": "H2A.R",
                    "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein", "nucleic": "Nucleic"}


        cxy=os.path.abspath(filex).split("round")[1].split("/")[1].split("/anal")[0]
        if cxy == "cenpa":
            sel1 = cenpa_chains[sel1x]
            sel2 = cenpa_chains[sel2x]
            print(sel1, sel2)
        elif cxy == "cenpa-h2az":
            sel1 = cenpa_h2az_chains[sel1x]
            sel2 = cenpa_h2az_chains[sel2x]
            print(sel1, sel2)
        elif cxy == "h3":
            sel1 = h3_chains[sel1x]
            sel2 = h3_chains[sel2x]
            print(sel1,sel2)

        path=os.getcwd()
        if not os.path.exists("00_graph/contact/heatmap"):
            os.makedirs("00_graph/contact/heatmap", exist_ok=True)
        pathx = path + "/00_graph/contact/heatmap"


        df= pd.read_csv(filex)
        pathh=os.path.abspath(filex)
        ro = pathh.split("MDx/")[1].split(f'/{cxy}')[0]
        file = os.path.basename(f'{filex}')
        chain = file.split(".csv")[0].split("chain")[1].split("_")[0]

        time=[]
        for i in range(1,len(df.columns)-1):
            time.append(df.columns[i].split("number_")[1])

        print(cxy,ro,chain)

        if chain == "J":
            df=df[::-1]
            df.reset_index(drop=True, inplace=True)
        if chain == "J" and cxy == "h3" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "I" and cxy == "h3" and ro == "round_2":
            df.resid = df.resid - 73

        if chain == "I" and cxy == "cenpa-h2az" and ro == "round_1":
            df.resid = df.resid - 73
        if chain == "J" and cxy == "cenpa-h2az" and ro == "round_1":
            df.resid = df.resid - 219


        if chain == "I" and  cxy == "cenpa" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "J" and  cxy == "cenpa" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "I" and  cxy == "cenpa-h2az" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "J" and  cxy == "cenpa-h2az" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "C":
            if cxy == "cenpa-h2az":
                df.resid = df.resid - 800
            #elif cx == "cenpa":
            #    aa.resid = aa.resid - 800



        df_ultimate=pd.DataFrame(columns=["resid", "number", "time"])
        for i in range(len(df)):
            lf_2 = [[df["resid"][i]] * len(time), df.iloc[i, 1:-1].tolist(), time]
            df_2=pd.DataFrame(lf_2).T
            df_2.columns = ["resid", "number", "time"]
            df_2.time = [int(i)/10 for i in df_2.time]
            df_ultimate = pd.concat([df_ultimate,df_2])




        df_ultimate.resid = df_ultimate.resid.astype(int)
        df_ultimate.time = df_ultimate.time.astype(int)
        df_ultimate.number = df_ultimate.number.astype(int)





        heatmap1_data = pd.pivot_table(df_ultimate, values='number',
                             index=['resid'],
                             columns='time')

        cmap = sns.cubehelix_palette(n_colors=80,dark=0,light=0.975,gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        #flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
        #cmap = ListedColormap(sns.color_palette(flatui).as_hex())

        #cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [ "#f9f5f7","#143943"])

        fig,ax=plt.subplots(1, 1,figsize=(25,25))
        #sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
        #            ax=ax)

        sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                    ax=ax, vmin=0, vmax=80)

        for _, spine in ax.spines.items():
            spine.set_visible(True)

        for t in ax.xaxis.get_ticklines(): t.set_color('silver')
        for t in ax.yaxis.get_ticklines(): t.set_color('silver')

        plt.title(f'Contact Number Per Residue - {sel1} to {sel2}',
                  fontdict=dict(weight='bold'),fontsize = 24)
        plt.ylabel("Residue Number", fontdict=dict(weight='bold'),fontsize = 20)
        plt.xlabel("Time Scale (ns)", fontdict=dict(weight='bold'),fontsize = 20)
        yticks=ax.get_yticklabels()
        ax.set_yticklabels(yticks,fontdict=dict(weight='bold',fontsize=14),rotation=0)

        xticks_max=max(df_ultimate.time)
        #xticks_min=min(df_ultimate.time)
        #xticks = np.arange(0, xticks_max, 50)
        xticks=np.arange(0, 502, 50)
        ax.xaxis.set_ticks(xticks)
        ax.set_xticklabels(xticks,fontdict=dict(fontsize=14))
        file_name = f'{sel1}_to_{sel2}'
        plt.savefig(f'{pathx}/{file_name}.png',bbox_inches='tight',dpi=300)

        plt.close()

    def barplot(filex,selection1,selection2):

    #path="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_2/h3/analysisNuc/02_contact"
    #pathr="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/"
    #os.chdir(path)

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        path=os.getcwd()
        if not os.path.exists("00_graph/contact/barplot"):
            os.makedirs("00_graph/contact/barplot", exist_ok=True)
        pathx = path + "/00_graph/contact/barplot"


        df= pd.read_csv(filex)

        time=[]
        for i in range(1,len(df.columns)-1):
            time.append(df.columns[i].split("number_")[1])


        ind = []
        for j in range(len(df)):
            if all(df.iloc[j, 1:] == 0):
                ind.append(j)

        df.drop(df.index[ind], inplace=True)

        df=df.reset_index(drop=True)


        df_ultimate=pd.DataFrame(columns=["resid", "number", "time"])
        for i in range(len(df)):
            lf_2 = [[df["resid"][i]] * len(time), df.iloc[i, 1:-1].tolist(), time]
            df_2=pd.DataFrame(lf_2).T
            df_2.columns = ["resid", "number", "time"]
            df_ultimate = pd.concat([df_ultimate,df_2])




        df_ultimate.resid = df_ultimate.resid.astype(int)
        df_ultimate.time = df_ultimate.time.astype(int)
        df_ultimate.number = df_ultimate.number.astype(int)



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

class contactperres():
    def multiple_heatmaps(pathq, selection1, selection2):
        pathv = os.getcwd()
        #pathb = pathv + "/" + pathq
        pathb = pathv + "/" + pathq.split("/")[1]

        roundsl = {"round_1": "601", "round_2": "Alpha Satellite", "round_3": "601L"}
        cenpa_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L",
                        "chainE": "CENPA.R",
                        "chainF": "H4.R", "chainG": "H2A.R",
                        "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                        "nucleic": "Nucleic"}
        cenpa_h2az_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.Z.L", "chainD": "H2B.L",
                             "chainE": "CENPA.R",
                             "chainF": "H4.R", "chainG": "H2A.Z.R",
                             "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                             "nucleic": "Nucleic"}
        h3_chains = {"chainA": "H3.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L", "chainE": "H3.R",
                     "chainF": "H4.R", "chainG": "H2A.R",
                     "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                     "nucleic": "Nucleic"}

        sel1x = selection1.replace(" ", "")
        sel2x = selection2.replace(" ", "")

        file_name = f'{sel1x}_to_{sel2x}'

        if not os.path.exists(f'{pathv}/00_graph/contactperres/multiple_heatmaps'):
            os.makedirs(f'{pathv}/00_graph/contactperres/multiple_heatmaps', exist_ok=True)
        pathx = f'{pathv}/00_graph/contactperres/multiple_heatmaps'

        #round = pathb.split("MDx/")[1]
        #round = pathq
        round = pathq.split("/")[1]
        cx = ["h3", "cenpa", "cenpa-h2az"]
        cx_rmsf = OrderedDict()
        sells = OrderedDict()

        for cxy in cx:

            filex = f'{pathv}/{round}/{cxy}/analysisNuc/06_contact_perres_merged/{file_name}/{file_name}_perres_num.csv'

            if cxy == "cenpa":
                sel1 = cenpa_chains[sel1x]
                sel2 = cenpa_chains[sel2x]
                print(sel1, sel2)
            elif cxy == "cenpa-h2az":
                sel1 = cenpa_h2az_chains[sel1x]
                sel2 = cenpa_h2az_chains[sel2x]
                print(sel1, sel2)
            elif cxy == "h3":
                sel1 = h3_chains[sel1x]
                sel2 = h3_chains[sel2x]
                print(sel1, sel2)

            sells[cxy] = [sel1, sel2]

            df = pd.read_csv(filex)
            pathh = os.path.abspath(filex)
            file = os.path.basename(f'{filex}')
            chain1 = file.split(".csv")[0].split("chain")[1].split("_")[0]
            chain2 = file.split(".csv")[0].split("_to_chain")[1].split("_")[0]

            df_rm = df.iloc[:, 1:]
            df_rm = df_rm.fillna(0)

            dfk = pd.DataFrame(df_rm.mean(axis=0))

            dfk = dfk.T

            r = re.compile("([a-zA-Z]+)([0-9]+)")


            donor_l = []
            for i, columns in enumerate(dfk.columns):  # add if statement
                donor_m = r.match(dfk.columns[i].split("_")[0])
                # donor_aa = donor_m.group(1)
                donor_num = donor_m.group(2)
                donor_num = int(donor_num)
                if chain1 == "C":
                    if cxy == "cenpa-h2az":
                        donor_num = donor_num - 800
                donor_l.append(donor_num)
            # donor_chainID = dfk.columns[0].split("_")[1].split(":")[0]

            acceptor_l = []
            for i, columns in enumerate(dfk.columns):  # add if statement
                acceptor_m = r.match(dfk.columns[i].split(":")[1].split("number_")[0])
                # acceptor_aa = acceptor_m.group(1)
                acceptor_num = acceptor_m.group(2)
                acceptor_num = int(acceptor_num)
                if chain2 == "C":
                    if cxy == "cenpa-h2az":
                        acceptor_num = acceptor_num - 800
                acceptor_l.append(acceptor_num)
            # acceptor_chainID = dfk.columns[0].split(":")[1].split("_")[1]

            donor_ls = sorted(set(donor_l))
            acceptor_ls = sorted(set(acceptor_l))

            fd = pd.DataFrame(columns=donor_ls, index=acceptor_ls)
            fd = fd.fillna(0)

            for i, columns in enumerate(dfk.columns):
                acceptor_m = r.match(dfk.columns[i].split(":")[1].split("_")[0])
                acceptor_num = acceptor_m.group(2)
                acceptor_num = int(acceptor_num)
                if chain2 == "C":
                    if cxy == "cenpa-h2az":
                        acceptor_num = acceptor_num - 800
                donor_m = r.match(dfk.columns[i].split("_")[0])
                donor_num = donor_m.group(2)
                donor_num = int(donor_num)
                if chain1 == "C":
                    if cxy == "cenpa-h2az":
                        donor_num = donor_num - 800
                fd.loc[acceptor_num, donor_num] = dfk.iloc[0, i]

            heatmap1_data = fd

            cx_rmsf[cxy] = heatmap1_data

        cmap = sns.cubehelix_palette(n_colors=10, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        fig, ax = plt.subplots(1, 3, figsize=(60, 25))

        for i, cxy in enumerate(cx):

            if i < 2:
                sns.heatmap(cx_rmsf[cxy], cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                            ax=ax[i], vmin=0, vmax=30, cbar=False)
            else:
                sns.heatmap(cx_rmsf[cxy], cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                            ax=ax[i], vmin=0, vmax=30)

            for _, spine in ax[i].spines.items():
                spine.set_visible(True)

            for t in ax[i].xaxis.get_ticklines(): t.set_color('silver')
            for t in ax[i].yaxis.get_ticklines(): t.set_color('silver')

            ax[i].set_title(f'{cxy.upper()} Complex - ({ sells[cxy][0]} to {sells[cxy][1]})',
                      fontdict=dict(weight='bold'), fontsize=24)
            ax[i].set_ylabel(f'Residue Number - {sells[cxy][1]}', fontdict=dict(weight='bold'), fontsize=20)
            ax[i].set_xlabel(f'Residue Number - {sells[cxy][0]}', fontdict=dict(weight='bold'), fontsize=20)
            yticks = ax[i].get_yticklabels()
            ax[i].set_yticklabels(yticks, fontdict=dict(weight='bold', fontsize=14), rotation=0)

            xticks = ax[i].get_xticklabels()
            ax[i].set_xticklabels(xticks, fontdict=dict(weight='bold', fontsize=14), rotation=0)
            plt.suptitle(f' Contact number at {roundsl[round]} Sequence', y=1.05, fontweight="bold", fontsize=30,
                         fontdict=dict(weight='bold'))

        fig.tight_layout()
        logging.info(f'writing to {pathx} ...')
        plt.savefig(f'{pathx}/{round}_{file_name}.png', bbox_inches='tight', dpi=300)
        plt.close()

    def heatmaps(filex, selection1, selection2):

        # path="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_2/h3/analysisNuc/02_contact"
        # pathr="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/"
        # os.chdir(path)

        sel1x = selection1.replace(" ", "")
        sel2x = selection2.replace(" ", "")

        cenpa_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L",
                        "chainE": "CENPA.R",
                        "chainF": "H4.R", "chainG": "H2A.R",
                        "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                        "nucleic": "Nucleic"}

        cenpa_h2az_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.Z.L", "chainD": "H2B.L",
                             "chainE": "CENPA.R",
                             "chainF": "H4.R", "chainG": "H2A.Z.R",
                             "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                             "nucleic": "Nucleic"}

        h3_chains = {"chainA": "H3.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L", "chainE": "H3.R",
                     "chainF": "H4.R", "chainG": "H2A.R",
                     "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                     "nucleic": "Nucleic"}

        cxy = os.path.abspath(filex).split("round")[1].split("/")[1].split("/anal")[0]
        if cxy == "cenpa":
            sel1 = cenpa_chains[sel1x]
            sel2 = cenpa_chains[sel2x]
            print(sel1, sel2)
        elif cxy == "cenpa-h2az":
            sel1 = cenpa_h2az_chains[sel1x]
            sel2 = cenpa_h2az_chains[sel2x]
            print(sel1, sel2)
        elif cxy == "h3":
            sel1 = h3_chains[sel1x]
            sel2 = h3_chains[sel2x]
            print(sel1, sel2)

        path = os.getcwd()
        if not os.path.exists("00_graph/contact_perres/heatmap"):
            os.makedirs("00_graph/contact_perres/heatmap", exist_ok=True)
        pathx = path + "/00_graph/contact_perres/heatmap"

        df = pd.read_csv(filex)
        pathh = os.path.abspath(filex)
        ro = pathh.split("MDx/")[1].split(f'/{cxy}')[0]
        file = os.path.basename(f'{filex}')
        chain1 = file.split(".csv")[0].split("chain")[1].split("_")[0]
        chain2 = file.split(".csv")[0].split("_to_chain")[1].split("_")[0]

        df_rm = df.iloc[:, 1:]
        df_rm = df_rm.fillna(0)

        dfk = pd.DataFrame(df_rm.mean(axis=0))

        dfk = dfk.T

        r = re.compile("([a-zA-Z]+)([0-9]+)")

        donor_l = []
        for i, columns in enumerate(dfk.columns):  # add if statement
            donor_m = r.match(dfk.columns[i].split("_")[0])
            # donor_aa = donor_m.group(1)
            donor_num = donor_m.group(2)
            donor_num = int(donor_num)
            if chain1 == "C":
                if cxy == "cenpa-h2az":
                    donor_num = donor_num - 800
            donor_l.append(donor_num)
        # donor_chainID = dfk.columns[0].split("_")[1].split(":")[0]

        acceptor_l = []
        for i, columns in enumerate(dfk.columns):  # add if statement
            acceptor_m = r.match(dfk.columns[i].split(":")[1].split("_")[0])
            # acceptor_aa = acceptor_m.group(1)
            acceptor_num = acceptor_m.group(2)
            acceptor_num = int(acceptor_num)
            if chain2 == "C":
                if cxy == "cenpa-h2az":
                    acceptor_num = acceptor_num - 800
            acceptor_l.append(acceptor_num)
        # acceptor_chainID = dfk.columns[0].split(":")[1].split("_")[1]

        donor_ls = sorted(set(donor_l))
        acceptor_ls = sorted(set(acceptor_l))

        fd = pd.DataFrame(columns=donor_ls, index=acceptor_ls)
        fd = fd.fillna(0)

        for i, columns in enumerate(dfk.columns):
            acceptor_m = r.match(dfk.columns[i].split(":")[1].split("_")[0])
            acceptor_num = acceptor_m.group(2)
            acceptor_num = int(acceptor_num)
            if chain2 == "C":
                if cxy == "cenpa-h2az":
                    acceptor_num = acceptor_num - 800
            donor_m = r.match(dfk.columns[i].split("_")[0])
            donor_num = donor_m.group(2)
            donor_num = int(donor_num)
            if chain1 == "C":
                if cxy == "cenpa-h2az":
                    donor_num = donor_num - 800
            fd.loc[acceptor_num, donor_num] = dfk.iloc[0, i]

        heatmap1_data = fd

        cmap = sns.cubehelix_palette(n_colors=10, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        fig, ax = plt.subplots(1, 1, figsize=(25, 25))

        sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                    ax=ax, vmin=0, vmax=10)

        for _, spine in ax.spines.items():
            spine.set_visible(True)

        for t in ax.xaxis.get_ticklines(): t.set_color('silver')
        for t in ax.yaxis.get_ticklines(): t.set_color('silver')

        plt.title(f'Contact Number Per Residue - {sel1} to {sel2}',
                  fontdict=dict(weight='bold'), fontsize=24)
        plt.ylabel(f'Residue Number - {sel1}', fontdict=dict(weight='bold'), fontsize=20)
        plt.xlabel(f'Residue Number - {sel2}', fontdict=dict(weight='bold'), fontsize=20)
        yticks = ax.get_yticklabels()
        ax.set_yticklabels(yticks, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        xticks = ax.get_xticklabels()
        ax.set_xticklabels(xticks, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        file_name = f'{sel1}_to_{sel2}'

        plt.savefig(f'{pathx}/{file_name}.png', bbox_inches='tight', dpi=300)

        plt.close()






class sasa():
    def sasa_heatmaps(filex, selection1, selection2):
        # path="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_2/h3/analysisNuc/02_contact"
        # pathr="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/"
        # os.chdir(path)

        sel1x = selection1.replace(" ", "")
        sel2x = selection2.replace(" ", "")

        cenpa_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L",
                        "chainE": "CENPA.R",
                        "chainF": "H4.R", "chainG": "H2A.R",
                        "chainH": "H2B.R", "chainI": "DNA I", "I": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                        "nucleic": "Nucleic", "all": ".", "chainIJ": "bp"}

        cenpa_h2az_chains = {"chainA": "CENPA.L", "chainB": "H4.L", "chainC": "H2A.Z.L", "chainD": "H2B.L",
                             "chainE": "CENPA.R",
                             "chainF": "H4.R", "chainG": "H2A.Z.R", "J": "DNA J",
                             "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                             "nucleic": "Nucleic", "all": ".", "chainIJ": "bp"}

        h3_chains = {"chainA": "H3.L", "chainB": "H4.L", "chainC": "H2A.L", "chainD": "H2B.L", "chainE": "H3.R",
                     "chainF": "H4.R", "chainG": "H2A.R", "I": "DNA I", "J": "DNA J",
                     "chainH": "H2B.R", "chainI": "DNA I", "chainJ": "DNA J", "protein": "Protein",
                     "nucleic": "Nucleic", "all": ".", "chainIJ": "bp"}

        cxy = os.path.abspath(filex).split("round")[1].split("/")[1].split("/anal")[0]
        if cxy == "cenpa":
            sel1 = cenpa_chains[sel1x]
            sel2 = cenpa_chains[sel2x]
            print(sel1, sel2)
        elif cxy == "cenpa-h2az":
            sel1 = cenpa_h2az_chains[sel1x]
            sel2 = cenpa_h2az_chains[sel2x]
            print(sel1, sel2)
        elif cxy == "h3":
            sel1 = h3_chains[sel1x]
            sel2 = h3_chains[sel2x]
            print(sel1, sel2)

        path = os.getcwd()
        if not os.path.exists("00_graph/sasa/heatmap"):
            os.makedirs("00_graph/sasa/heatmap", exist_ok=True)
        pathx = path + "/00_graph/sasa/heatmap"

        df = pd.read_csv(filex)
        df.drop('mean', axis=1, inplace=True)

        pathh = os.path.abspath(filex)
        ro = pathh.split("MDx/")[1].split(f'/{cxy}')[0]
        file = os.path.basename(f'{filex}')
        chain = file.split(".csv")[0].split("chain")[1].split("_")[0]

        time = []
        for i in range(1, len(df.columns)):
            time.append(df.columns[i].split("sasa_")[1])

        print(cxy, ro, chain)

        if chain == "J":
            df = df[::-1]
            df.reset_index(drop=True, inplace=True)
        if chain == "J" and cxy == "h3" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "I" and cxy == "h3" and ro == "round_2":
            df.resid = df.resid - 73

        if chain == "I" and cxy == "cenpa-h2az" and ro == "round_1":
            df.resid = df.resid - 73
        if chain == "J" and cxy == "cenpa-h2az" and ro == "round_1":
            df.resid = df.resid - 219

        if chain == "I" and cxy == "cenpa" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "J" and cxy == "cenpa" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "I" and cxy == "cenpa-h2az" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "J" and cxy == "cenpa-h2az" and ro == "round_2":
            df.resid = df.resid - 73
        if chain == "C":
            if cxy == "cenpa-h2az":
                df.resid = df.resid - 800

        # heatmap1_data = pd.pivot_table(df, values='number',
        #                               index=['resid'],
        #                               columns='time')

        heatmap1_data = df

        if filex.endswith("_notbp.csv"):
            heatmap1_data.set_index(df.loc[:, "resid"], inplace=True)
            heatmap1_data.drop('resid', axis=1, inplace=True)

        if filex.endswith("_bp.csv"):
            heatmap1_data.set_index(df.loc[:, "resid_bp"], inplace=True)
            heatmap1_data.drop('resid_bp', axis=1, inplace=True)
        heatmap1_data.columns = time

        cmap = sns.cubehelix_palette(n_colors=400, dark=0, light=0.975, gamma=1.4,
                                     rot=-0.2855, as_cmap=True)

        # flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
        # cmap = ListedColormap(sns.color_palette(flatui).as_hex())

        # cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [ "#f9f5f7","#143943"])

        fig, ax = plt.subplots(1, 1, figsize=(25, 25))
        # sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
        #            ax=ax)

        sns.heatmap(heatmap1_data, cmap=cmap, yticklabels=True, cbar_kws={"shrink": 0.5},
                    ax=ax, vmin=0, vmax=400)

        for _, spine in ax.spines.items():
            spine.set_visible(True)

        for t in ax.xaxis.get_ticklines(): t.set_color('silver')
        for t in ax.yaxis.get_ticklines(): t.set_color('silver')

        plt.title(f'SASA - {sel2}',
                  fontdict=dict(weight='bold'), fontsize=24)
        plt.ylabel("Residue Number", fontdict=dict(weight='bold'), fontsize=20)
        plt.xlabel("Time Scale (ns)", fontdict=dict(weight='bold'), fontsize=20)
        yticks = ax.get_yticklabels()
        ax.set_yticklabels(yticks, fontdict=dict(weight='bold', fontsize=14), rotation=0)

        # xticks_max = max(time)
        # xticks_max = 502
        # xticks_min=min(df_ultimate.time)
        # xticks = np.arange(0, xticks_max, 50)
        xticks = np.arange(0, 5002, 500)
        ax.xaxis.set_ticks(xticks)
        xticks_la = np.arange(0, 502, 50)
        ax.set_xticklabels(xticks_la, fontdict=dict(fontsize=14))
        file_name = f'{sel2}_{ro}'
        # plt.show()
        plt.savefig(f'{pathx}/{file_name}.png', bbox_inches='tight', dpi=300)

        plt.close()


class rmsf():
    def barplot(round,title,system_type):

        types = ['system', 'protein', 'dna']
        if system_type not in types:
            raise ValueError("Invalid system type. Expected one of: %s" % types)

        path=os.getcwd()
        if not os.path.exists(f'{path}/{round}/00_graph/rmsfPerResidue/barplot'):
            os.makedirs(f'{path}/{round}/00_graph/rmsfPerResidue/barplot', exist_ok=True)
        pathx = f'{path}/{round}/00_graph/rmsfPerResidue/barplot'


        if system_type == "system":
            chains = ["A","B","C","D","E","F","G","H","I","J"]
        elif system_type == "protein":
            chains = ["A","B","C","D","E","F","G","H"]
        elif system_type == "dna":
            chains = ["I","J"]

        cx =["h3","cenpa","cenpa-h2az"]
        cx_rmsf = OrderedDict()
        for cxx in cx:
            me_rmsf = OrderedDict()
            for chain in chains:
                #print(cxx,chain)
                aa = pd.read_csv(
                    f'{path}/{round}/{cxx}/analysisNuc/03_rmsfPerResidue/rmsf_chain{chain}.csv')
                me_rmsf[chain] = aa
            cx_rmsf[cxx] = me_rmsf



        fig, ax_box = plt.subplots(len(chains),3, figsize=(12,len(chains)*3), sharey=True)

        cxl =["H3","CENPA","CENPA-H2AZ"]
        plot_ = OrderedDict()
        #plot_x = OrderedDict()
        colors = ["#ccd9b4","#6F81A6","#D9B4CC"]
        for j, cxy in enumerate(cx):
            for i, chainy in enumerate(chains):
                plot_[cxy,chainy]=sns.barplot(x="resid", y="rmsf", data=cx_rmsf[cxy][chainy], color=colors[j],ax=ax_box[i,j])
                for ind, label in enumerate(plot_[cxy,chainy].get_xticklabels()):
                    if ind % 10 == 0:  # every 10th label is kept
                        label.set_visible(True)
                    else:
                        label.set_visible(False)
                ax_box[i,j].set_title(f'{chainy}', fontdict=dict(weight='bold'),fontsize = 16)
                ax_box[i,j].tick_params(axis='x', labelsize=18, rotation=90)
                ax_box[i, j].tick_params(axis='y', labelsize=20)
                ax_box[i,j].set_xlabel('')
                ax_box[i,0].set_ylabel('RMSF', fontdict=dict(weight='bold'),fontsize = 16)
                ax_box[i, 1].set_ylabel('')
                ax_box[i, 2].set_ylabel('')

        #plt.ylabel('Contact Number',fontsize = 12)
            ax_box[i,j].set_xlabel(f'{cxl[j]}',fontsize = 16, fontdict=dict(weight='bold'))
        #plt.legend(('vmd', 'python'), loc='upper right')
        plt.suptitle(f'RMSF per Residue - {title}',y=1.02,fontweight="bold", fontsize=24, fontdict=dict(weight='bold'))
        fig.tight_layout()
        plt.savefig(f'{pathx}/{system_type}_rmsfperres_{round}.png',bbox_inches='tight',dpi=300)
        plt.close()


