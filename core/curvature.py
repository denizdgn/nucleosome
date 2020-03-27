#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd
import pathos
from pathos.multiprocessing import ProcessingPool as Pool
import sys
from collections import OrderedDict
import logging
import io
import numpy as np
import glob
import numpy.linalg as la
import matplotlib.pylab as plt
import seaborn as sns
import itertools
import signal


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')



class curvatured():
    def final_call(target,tchainID1,tchainID2):
        logging.info(f'Starting to analyze {target} ...')
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        pathx = os.getcwd()

        def read_pdb(target):
            """
            :param target: single pdb file
            :return:
            """

            #path = os.getcwd()
            f = open(target, "r")
            x = f.read()
            # Read a table of fixed-width formatted lines into DataFrame
            colnames = ["type", "atom_num", "atom_name", "idk", "res_name", "chainID", "res_num", "idk_code", "X", "Y", "Z",
                        "occupancy",
                        "temp_fac", "segid", "element_sym", "charge"]


            widths = [6, 5, 5, 1, 3, 2, 4, 1, 11, 8, 8, 6, 6, 10, 2, 2]
            df = pd.read_fwf(io.StringIO(x), header=None, widths=widths, names=colnames)


            df_filt_1 = df[df.type == "ATOM"]
            df_filt_2 = df[df.type == "HETATM"]
            df_filt = df_filt_1.append(df_filt_2)
            df_filt.drop(["occupancy", "temp_fac", "segid", "element_sym", "charge"], axis=1, inplace=True)
            df_filt = df_filt.fillna('')

            return df_filt



        def get_residue(df_filt, tchainID,tresnum):
            chainID_filt = (df_filt.chainID == str(tchainID))
            res_num_filt = (df_filt.res_num == int(tresnum))
            residue = df_filt[chainID_filt & res_num_filt]
            return residue


        def get_base_center(residue):
            """
            To find geometric center of a residue
            :param residue: residue in pdb format
            :return: geometric center
            """
            resarray = np.asarray(residue[["X","Y","Z"]])

            com = np.zeros(3)
            natoms = len(resarray)

            for i,j in itertools.product(range(0, natoms),range(0, 3)):
                com[j] += resarray[i][j]

            com /= natoms
            return com


        def get_bp_center(df_filt, tchainID1, tchainID2, tresnum1, tresnum2):
            base1 = get_residue(df_filt, tchainID1,tresnum1)
            base2 = get_residue(df_filt, tchainID2, tresnum2)
            return 0.5 * (get_base_center(base1) + get_base_center(base2))



        def get_bp_centers(df_filt,tchainID1, tchainID2):
            A=df_filt[df_filt.chainID == tchainID1].res_num
            resA=list(OrderedDict.fromkeys(A).keys())


            B=df_filt[df_filt.chainID == tchainID2].res_num
            resB=list(OrderedDict.fromkeys(B).keys())
            resB.reverse()

            xl=[]
            for a,b in zip(resA,resB):
                xl.append(get_bp_center(df_filt, tchainID1, tchainID2, a, b))
            xlx = np.array(xl)

            return xlx


        def write_xyz(com, filename, color_for_element):
            """

            """
            outfile = open(filename, 'a')

            ncom = com.shape[0]

            outfile.write(str(ncom) + "\n")
            #time = str(target.split("md_")[1].split(".pd")[0])
            #outfile.write(str(time) + "\n")
            for i in range(0, ncom):

                outfile.write(color_for_element + " " + str(com[i][0]) + " " + str(com[i][1]) + " " + str(com[i][2]) + "\n")

            outfile.close()


        def get_bp_bp_distances(com):
            ndist = com.shape[0] - 2
            distances = np.zeros(ndist)

            for i in range(0, ndist):
                distances[i] = la.norm(com[i] - com[i + 1])

            return distances

        def py_ang(v1, v2):
            """ Returns the angle in radians between vectors 'v1' and 'v2'    """
            cosang = np.dot(v1, v2)
            sinang = la.norm(np.cross(v1, v2))

            return np.arctan2(sinang, cosang)


        def get_bp_bp_angles(com):
            nangles = com.shape[0] - 2
            angles = np.zeros(nangles)

            for i in range(0, nangles):
                j = i + 1
                angles[i] = py_ang(com[j + 1] - com[j], com[j - 1] - com[j])

                if angles[i] < np.pi / 2:
                    angles[i] = np.pi - angles[i]

            return np.rad2deg(angles)


        def get_radius_of_curvature(com):
            nradii = com.shape[0] - 2
            radii = np.zeros(nradii)

            for i in range(0, nradii):
                j = i + 1

                d1 = la.norm(com[j] - com[j + 1])
                d2 = la.norm(com[j] - com[j - 1])
                d3 = la.norm(com[j + 1] - com[j - 1])

                # theta = py_ang(com[i+1]-com[i], com[i-1]-com[i])

                # Find area by Heron's formula
                s = 0.5 * (d1 + d2 + d3)

                A = np.sqrt(s * (s - d1) * (s - d2) * (s - d3))

                radii[i] = (d1 * d2 * d3) / (4. * A)

            signal.signal(signal.SIGINT, handler)

            # return radii, 1. / radii
            return radii


        def smoothen_com(com, window):
            """

            Use averaging over a window of neighboring centers

            """

            nbp = com.shape[0]

            halfwindow = int((window - 1) / 2)

            scom = np.zeros([nbp - 2 * halfwindow, 3])

            for i in range(halfwindow, nbp - halfwindow):
                for j in range(-halfwindow, halfwindow + 1):
                    scom[i - halfwindow] += com[i + j]
                    signal.signal(signal.SIGINT, handler)

            return scom / window


        def get_bending_energy(distances, radii):
            energy = 0  # np.zeros(radii)

            for i in range(1, len(distances) - 1):
                arc = 0.5 * (distances[i] + distances[i - 1])
                signal.signal(signal.SIGINT, handler)

                energy += 0.5 * 500 * arc / radii[i] ** 2

            return energy



        time = str(target.split("md_")[1].split(".pd")[0])


        df_filt = read_pdb(target)



        #pathx = '/Users/denizdogan/Desktop/xx'
        pathx = os.getcwd()
        if not os.path.exists(pathx+"/07_nuc_com_xyz"):
            os.makedirs(pathx+"/07_nuc_com_xyz")
        pathxz = pathx + "/07_nuc_com_xyz"


        nuc_com_sm1=get_bp_centers(df_filt,tchainID1, tchainID2)
        write_xyz(nuc_com_sm1, f'{pathxz}/{time}t_nuc_com_sm1.xyz', 'O')

        nuc_com_sm3 = smoothen_com(nuc_com_sm1, 3)
        write_xyz(nuc_com_sm3, f'{pathxz}/{time}t_nuc_com_sm3.xyz', 'P')

        nuc_com_sm5 = smoothen_com(nuc_com_sm1, 5)
        write_xyz(nuc_com_sm5, f'{pathxz}/{time}t_nuc_com_sm5.xyz', 'S')

        nuc_com_sm7 = smoothen_com(nuc_com_sm1, 7)
        write_xyz(nuc_com_sm7, f'{pathxz}/{time}t_nuc_com_sm7.xyz', 'N')

        nuc_com_sm9 = smoothen_com(nuc_com_sm1, 9)
        write_xyz(nuc_com_sm9, f'{pathxz}/{time}t_nuc_com_sm9.xyz', 'C')


        nuc_rc_sm1 = get_radius_of_curvature(nuc_com_sm1)
        nuc_rc_sm3 = get_radius_of_curvature(nuc_com_sm3)
        nuc_rc_sm5 = get_radius_of_curvature(nuc_com_sm5)
        nuc_rc_sm7 = get_radius_of_curvature(nuc_com_sm7)
        nuc_rc_sm9 = get_radius_of_curvature(nuc_com_sm9)



        #os.chdir("..")
        #pathx='/Users/denizdogan/Desktop/xx'
        pathx = os.getcwd()
        if not os.path.exists(pathx+"/08_nuc_rc_sm"):
            os.makedirs(pathx+"/08_nuc_rc_sm")
        pathxx = pathx + "/08_nuc_rc_sm"


        j=1
        for value in [nuc_rc_sm1,nuc_rc_sm3,nuc_rc_sm5,nuc_rc_sm7,nuc_rc_sm9]:
            #print(value)
            df = pd.DataFrame(value).T
            df["time"] = time
            df["mean"] = np.mean(value)
            df.to_csv(f'{pathxx}/{time}t_nuc_rc_sm{j}.csv', index=None)
            j+=2
            signal.signal(signal.SIGINT, handler)


        nuc_d_sm1 = get_bp_bp_distances(nuc_com_sm1)
        nuc_d_sm3 = get_bp_bp_distances(nuc_com_sm3)
        nuc_d_sm5 = get_bp_bp_distances(nuc_com_sm5)
        nuc_d_sm7 = get_bp_bp_distances(nuc_com_sm7)
        nuc_d_sm9 = get_bp_bp_distances(nuc_com_sm9)


        nuc_e = np.zeros([5])

        nuc_e[0] = get_bending_energy(nuc_d_sm1, nuc_rc_sm1)
        nuc_e[1] = get_bending_energy(nuc_d_sm3, nuc_rc_sm3)
        nuc_e[2] = get_bending_energy(nuc_d_sm5, nuc_rc_sm5)
        nuc_e[3] = get_bending_energy(nuc_d_sm7, nuc_rc_sm7)
        nuc_e[4] = get_bending_energy(nuc_d_sm9, nuc_rc_sm9)


        #pathx = os.getcwd()
        #pathx = '/Users/denizdogan/Desktop/xx'
        if not os.path.exists(pathx+"/09_nuc_e"):
            os.makedirs(pathx+"/09_nuc_e")
        pathxy = pathx + "/09_nuc_e"

        df_nuc_e = pd.DataFrame(nuc_e).T
        df_nuc_e["time"] = time
        df_nuc_e.to_csv(f'{pathxy}/{time}t_nuc_e.csv', index=None)


    def combine_xyz():
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")
        windows = np.linspace(1, 9, 5, dtype=int)
        pathx = os.getcwd()

        if not os.path.exists(pathx+"/07_nuc_com_xyz_merged"):
            os.makedirs(pathx+"/07_nuc_com_xyz_merged")
        pathxz = pathx + "/07_nuc_com_xyz_merged"

        for i in windows:
            df3 = pd.DataFrame()
            logging.info(f'Merging *nuc_com_xyz{i}.xyz files ...')
            files_xyz = glob.glob(pathx + f'/07_nuc_com_xyz/*{i}.xyz')
            for j in range(len(files_xyz)):
                file=pathx + f'/07_nuc_com_xyz/{j}t_nuc_com_sm{i}.xyz'
                df = pd.read_csv(file, header=None)
                df2 = pd.concat([df.iloc[:1], pd.DataFrame([0]), df.iloc[1:]]).reset_index(drop=True)
                df3=df3.append(df2)
            signal.signal(signal.SIGINT, handler)
            df3.to_csv(pathxz + f'/nuc_com_xyz{i}.xyz', index=None, header=None)



    def combine_rc():
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")
        windows = np.linspace(1, 9, 5, dtype=int)
        pathx = os.getcwd()

        if not os.path.exists(pathx+"/08_nuc_rc_sm_merged"):
            os.makedirs(pathx+"/08_nuc_rc_sm_merged")
        pathxz = pathx + "/08_nuc_rc_sm_merged"
        for i in windows:
            df3 = pd.DataFrame()
            logging.info(f'Merging *nuc_com_xyz{i}.xyz files ...')
            files_csv = glob.glob(pathx + f'/08_nuc_rc_sm/*{i}.csv')
            for j in range(len(files_csv)):
                file=pathx + f'/08_nuc_rc_sm/{j}t_nuc_rc_sm{i}.csv'
                df = pd.read_csv(file)
                df3=df3.append(df)
            signal.signal(signal.SIGINT, handler)
            df3.to_csv(pathxz + f'/nuc_rc_sm{i}.csv', index=None)

    def combine_nuce():

        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        pathx = os.getcwd()
        if not os.path.exists(pathx+"/09_nuc_e_merged"):
            os.makedirs(pathx+"/09_nuc_e_merged")
        pathxz = pathx + "/09_nuc_e_merged"

        df3 = pd.DataFrame()
        logging.info(f'Merging *t_nuc_e.csv files ...')
        files_csv = glob.glob(pathx + f'/09_nuc_e/*.csv')
        for file in files_csv:
            df = pd.read_csv(file)
            df3=df3.append(df)
            signal.signal(signal.SIGINT, handler)
        dff=df3.sort_values(by=['time'])
        dff.to_csv(pathxz + f'/nuc_e.csv', index=None)


