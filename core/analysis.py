import os
import pandas as pd
import vmd
import sys
from vmd import *
from collections import OrderedDict
import logging
import signal
import pathos
import glob
import shutil
import itertools
import io


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')

def paste(x: 'donor_resnm', y: 'donor_resid', a=None, b=None,sep="") :
    """

    :param x: donor_resnm
    :param y: donor_resid
    :param a: acceptor_resnm
    :param b: acceptor_resid
    :return: pasted variables as string
    Ex: paste("a","b") // out: 'ab'
    """
    if a != None :
        return str (x) + str (y) + str(sep) + str (a) + str (b)
    else :
        return str (x)+ str(sep) + str (y)





class contact():
    def combine_contact(selection1,selection2):
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        pathx = os.getcwd()
        pathx_toreadz = pathx + f'/06_contact/{sel1}_to_{sel2}'

        folder = f'{sel1}_to_{sel2}'
        logging.info(f'Running in folder: {folder} ...')
        #pathx_toreadz = pathx_toread + f'/{folder}'
        if not os.path.exists(pathx + f'/06_contact_merged/{folder}'):
            os.makedirs(pathx + f'/06_contact_merged/{folder}', exist_ok=True)
        pathx_towrite = pathx + f'/06_contact_merged/{folder}'

        df3 = pd.DataFrame()
        logging.info(f'Merging files ...')
        #get all files in a folder
        files_csv = glob.glob(pathx_toreadz + f'/*.csv')
        for time in range(len(files_csv)):
            df = pd.read_csv(f'{pathx_toreadz}/{sel1}_to_{sel2}_t{time}.csv')
            df3[df.columns[1]] = df.iloc[:,1]
            #signal.signal(signal.SIGINT, handler)
        df3["mean"] = df3.iloc[:, 1:].mean(axis=1)
        df4=pd.concat([pd.DataFrame(df.iloc[:,0]),df3],axis=1)
        df4.to_csv(f'{pathx_towrite}/{folder}.csv', index=None)

    def combine_contact_perres(selection1, selection2):
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        signal.signal(signal.SIGINT, handler)

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        pathx = os.getcwd()
        pathx_toreadz = pathx + f'/06_contact_perres/{sel1}_to_{sel2}'

        folder = f'{sel1}_to_{sel2}'
        logging.info(f'Running in folder: {folder} ...')
        # pathx_toreadz = pathx_toread + f'/{folder}'
        if not os.path.exists(pathx + f'/06_contact_perres_merged/{folder}'):
            os.makedirs(pathx + f'/06_contact_perres_merged/{folder}', exist_ok=True)
        pathx_towrite = pathx + f'/06_contact_perres_merged/{folder}'

        df3 = pd.DataFrame()
        logging.info(f'Merging files ...')
        # get all files in a folder
        files_csv = glob.glob(pathx_toreadz + f'/*perres_t*.csv')

        for time in range(0, 5003):
            try:
                df = pd.read_csv(f'{pathx_toreadz}/{sel1}_to_{sel2}_perres_t{time}.csv')

            except FileNotFoundError:
                df = pd.DataFrame(columns=['donor_chain', 'acceptor_chain', 'donor_resnm', 'acceptor_resnm',
                                           'donor_resid', 'acceptor_resid', 'donor_atom', 'acceptor_atom', 'donor',
                                           'donorC', 'acceptor', 'acceptorC', 'donor_acceptor', 'specificity',
                                           'time'])
                df.loc[0, "time"] = time

            df3 = df3.append(df)

        df3.to_csv(f'{pathx_towrite}/{folder}_perres.csv', index=None)

        df4 = pd.DataFrame()
        logging.info(f'Merging files ...')
        # get all files in a folder
        files_csv = glob.glob(pathx_toreadz + f'/*perres_num_t*.csv')
        for time in range(0, 5003):
            try:
                df2 = pd.read_csv(f'{pathx_toreadz}/{sel1}_to_{sel2}_perres_num_t{time}.csv')
                df2.insert(0, 'time', time)

            except FileNotFoundError:
                df2 = pd.DataFrame(columns=["time"])
                df2.loc[0, "time"] = time

            df4 = df4.append(df2)

        df4.to_csv(f'{pathx_towrite}/{folder}_perres_num.csv', index=None)


    def contact(pdb,selection1,selection2,cutoff):
        """
        usage: nucleus contact -p npt.pdb -s1 "chain A" -s2 "nucleic"
        :param pdb: Give traj file as .pdb
        :param selection1: First group to find contacts, i.e. "chain A", "nucleic"
        :param selection2: Second group to find contacts, i.e. "chain A", "nucleic"
        :return:
        """
        def handler():
            raise KeyboardInterrupt("CTRL-C!")

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")


        path=os.getcwd()
        os.chdir(path)
        if not os.path.exists(f'06_contact/{sel1}_to_{sel2}'):
            os.makedirs(f'06_contact/{sel1}_to_{sel2}', exist_ok=True)
        pathx = path + f'/06_contact/{sel1}_to_{sel2}'

        logging.info(f'Finding contacts from {sel1} to {sel2} ...')

        time = str(pdb.split("md_")[1].split(".pd")[0])
        mol = vmd.molecule.load('pdb', pdb)
        df_res_cont_num = pd.DataFrame()
        signal.signal(signal.SIGINT, handler)
        logging.info(f'Running on {pdb}')

        histone_chain = atomsel(f'{selection1} and noh')
        histone_chain.update()

        dna_chain = atomsel(f'{selection2} and noh')
        dna_chain.update()



        resids=list(OrderedDict.fromkeys(histone_chain.resid).keys())
        df_res_cont_num["resid"]=resids
        #iterate over residues
        contact_number = []
        for j in range(len(resids)):
            reside_atoms = atomsel(f'{selection1} and resid \"{resids[j]}\" and noh')
            reside_atoms.update()
            a=reside_atoms.contacts(selection=dna_chain,cutoff=float(cutoff))
            contact_number.append(len(a[1]))
        df_res_cont_num[f'contact_number_{time}'] = contact_number
        df_res_cont_num.to_csv(f'{pathx}/{sel1}_to_{sel2}_t{time}.csv',
                                   index=None)

    def contactperres(pdb, selection1, selection2, cutoff):
        f = open(f'{pdb}', "r")
        x = f.read()
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
        df_filt.atom_num = df_filt.atom_num.astype(int)
        time = str(pdb.split("md_")[1].split(".pd")[0])

        mol = vmd.molecule.load('pdb', f'{pdb}')

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        path = os.getcwd()
        os.chdir(path)
        if not os.path.exists(f'06_contact_perres/{sel1}_to_{sel2}'):
            os.makedirs(f'06_contact_perres/{sel1}_to_{sel2}', exist_ok=True)
        pathx = path + f'/06_contact_perres/{sel1}_to_{sel2}'
        logging.info(f'Finding contacts from {sel1} to {sel2} ...')

        select1 = atomsel(f'{selection1}', frame=1)
        resids_select1 = list(OrderedDict.fromkeys(select1.resid).keys())

        column_names = ['donor_chain', 'acceptor_chain', 'donor_resnm', 'acceptor_resnm', 'donor_resid',
                        'acceptor_resid', 'donor_atom', 'acceptor_atom']
        dfx = pd.DataFrame(columns=column_names)
        for resid1 in resids_select1:
            reside_atoms = atomsel(f'{selection1} and resid \"{resid1}\" and noh')
            reside_atoms.update()
            dna_resid = atomsel(f'{selection2} and noh')
            dna_resid.update()
            a = reside_atoms.contacts(selection=dna_resid, cutoff=float(cutoff))

            for i in range(len(a[0])):
                dict_l = {"donor_chain": df_filt.iloc[a[0][i], :].chainID,
                          "acceptor_chain": df_filt.iloc[a[1][i], :].chainID,
                          "donor_resnm": df_filt.iloc[a[0][i], :].res_name,
                          "acceptor_resnm": df_filt.iloc[a[1][i], :].res_name,
                          "donor_resid": df_filt.iloc[a[0][i], :].res_num,
                          "acceptor_resid": df_filt.iloc[a[1][i], :].res_num,
                          "donor_atom": df_filt.iloc[a[0][i], :].atom_name,
                          "acceptor_atom": df_filt.iloc[a[1][i], :].atom_name
                          }

                dfx = dfx.append(pd.DataFrame([dict_l]))
                del dict_l
            del a

        df_table = dfx
        df_table.donor_resid = df_table.donor_resid.astype(int)
        df_table.acceptor_resid = df_table.acceptor_resid.astype(int)



        donor_list=df_table.apply(lambda x: x['donor_resnm'] + str(x['donor_resid']), axis=1)
        df_table.loc[:, 'donor'] = donor_list

        donorC_list=df_table.apply(lambda x,sep="_": x['donor'] +"_"+ str(x['donor_chain']), axis=1)
        df_table.loc[:, 'donorC'] = donorC_list



        acceptor_list=df_table.apply(lambda x: x['acceptor_resnm'] + str(x['acceptor_resid']), axis=1)
        df_table.loc[:, 'acceptor'] = acceptor_list

        acceptorC_list=df_table.apply(lambda x,sep="_": x['acceptor'] +"_"+ str(x['acceptor_chain']), axis=1)
        df_table.loc[:, 'acceptorC'] = acceptorC_list


        donor_acceptor_list=df_table.apply(lambda x,sep=":": x['donorC'] +sep+ str(x['acceptorC']), axis=1)
        df_table.loc[:, 'donor_acceptor'] = donor_acceptor_list


        df_table.index = range(len(df_table))

        # ADD SPECIFICITY
        non_spp_atoms = ["O2P", "O1P", "N", "O", "OC1", "OC2", "O4'", "O5'", "O3'", "H", "HA",
                         "P", "C5'", "C4'", "C3'", "C2'", "C1'","CA","C","CB"]
        x = []
        for i in range(len(df_table)):
            if (df_table["acceptor_atom"][i] or df_table["donor_atom"][i]) in non_spp_atoms:
                x.append("non-specific")
            else:
                x.append("specific")
        df_table.loc[:, 'specificity'] = x
        df_table.loc[:, 'time'] = time

        df_table.to_csv(f'{pathx}/{sel1}_to_{sel2}_perres_t{time}.csv', index=None)

        xx = list(OrderedDict.fromkeys(df_table.donor_acceptor).keys())

        org = OrderedDict()
        for i, key in enumerate(xx):
            org[key] = sum(df_table.donor_acceptor == xx[i])


        org_df = pd.DataFrame([org])
        org_df.to_csv(f'{pathx}/{sel1}_to_{sel2}_perres_num_t{time}.csv', index=None)


def rmsfperresid(pdb,selection1,first,last,step):
    """
    :param pdb: Give traj file as .pdb
    :param selection1: To select a group of atoms, i.e "chain A"
    :param first: First frame to use, default is first frame
    :param last: Last frame to use, default is last frame
    :param step: To select every Nth frame
    :return: Write to file
    """
    logging.info(f'Starting..')

    path=os.getcwd()
    os.chdir(path)

    if not os.path.exists("03_rmsfPerResidue"):
        os.makedirs("03_rmsfPerResidue")
    pathx = path + "/03_rmsfPerResidue"

    mol = vmd.molecule.load('pdb', pdb)
    chain = atomsel(f'{selection1} and noh')
    chain.update()

    resids = list(OrderedDict.fromkeys(chain.resid).keys())

    df_res_rmsf = pd.DataFrame()
    df_res_rmsf["resid"] = resids

    rmsf = chain.rmsfperresidue(first=first,last=last,step=step)
    df_res_rmsf["rmsf"] = rmsf
    sel1 = selection1.replace(" ", "")
    logging.info(f'Writing ..')
    df_res_rmsf.to_csv(pathx + "/rmsf_" + sel1 +".csv",
                           index=None)





class sasa():
    def combine_sasa(selection1,selection2,bp):
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        pathx = os.getcwd()
        pathx_toreadz = pathx + f'/04_sasa/{sel1}_to_{sel2}'

        folder = f'{sel1}_to_{sel2}'
        logging.info(f'Running in folder: {folder} ...')


        if not os.path.exists(pathx + f'/04_sasa_merged/{folder}'):
            os.makedirs(pathx + f'/04_sasa_merged/{folder}', exist_ok=True)
        pathx_towrite = pathx + f'/04_sasa_merged/{folder}'

        df3 = pd.DataFrame()
        logging.info(f'Merging sasa files ...')

        #get all files in a folder
        if (bp):
            files_csv = glob.glob(pathx_toreadz + f'/*bp.csv')
            for time in range(len(files_csv)):
                df = pd.read_csv(f'{pathx_toreadz}/{sel1}_{sel2}_t{time}_bp.csv')
                df3[df.columns[1]] = df.iloc[:,1]
                signal.signal(signal.SIGINT, handler)
            df = pd.read_csv(f'{pathx_toreadz}/{sel1}_{sel2}_t83_bp.csv')
            df3["mean"] = df3.iloc[:, 0:].mean(axis=1)
            df4=pd.concat([pd.DataFrame(df.iloc[:,0]),df3],axis=1)
            logging.info(f'writing to file {folder}')
            df4.to_csv(f'{pathx_towrite}/{folder}_bp.csv', index=None)
        else:
            files_csv = glob.glob(pathx_toreadz + f'/*notbp.csv')
            for time in range(len(files_csv)):
                df = pd.read_csv(f'{pathx_toreadz}/{sel1}_{sel2}_t{time}_notbp.csv')
                df3[df.columns[1]] = df.iloc[:,1]
                #signal.signal(signal.SIGINT, handler)
            df = pd.read_csv(f'{pathx_toreadz}/{sel1}_{sel2}_t83_notbp.csv')
            df3["mean"] = df3.iloc[:, 0:].mean(axis=1)
            df4=pd.concat([pd.DataFrame(df.iloc[:,0]),df3],axis=1)
            logging.info(f'writing to file {folder}')
            df4.to_csv(f'{pathx_towrite}/{folder}_notbp.csv', index=None)

    def sasa(pdb,selection1:str,selection2:str,srad:float,bp):
        """
        :param pdb: *pdb trajectory
        :param selection1: everything in  the system (protein+DNA)
        :param selection2: chain IDs
        :param srad: Solvent radius
        :param bp: boolean - base pair option
        :return:
        """
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        class Error(Exception):
            pass
        class basepairError(Error):
            pass
        class nonNucleotideError(Error):
            pass

        path=os.getcwd()

        sel1 = selection1.replace(" ", "")
        sel2 = selection2.replace(" ", "")

        if not os.path.exists(f'04_sasa/{sel1}_to_{sel2}'):
            os.makedirs(f'04_sasa/{sel1}_to_{sel2}', exist_ok=True)
        pathx = path + f'/04_sasa/{sel1}_to_{sel2}'

        logging.info(f'Reading {pdb} .. ')
        mol = vmd.molecule.load('pdb', pdb)
        time = str(pdb.split("md_")[1].split(".pd")[0])

        select1 = atomsel(f'{selection1}')
        select2 = atomsel(f'{selection2}')

        try:
            if (len(set(select2.chain)) != 1):
                    if (all(select2.nucleic)):
                        if (bp):
                            pass
                        else:
                            raise basepairError
                    else:
                        raise nonNucleotideError

        except basepairError:
            print("\n\n!!!! basepairError: Use -bp option to calculate SASA for nucleotide base pairs !!!!")
            sys.exit(1)

        except nonNucleotideError:
            print("\n\n!!!! nonNucleotideError: base pairs - SASA calculations only enabled for nucleotides !!!!")
            sys.exit(1)


        resids_select2 = list(OrderedDict.fromkeys(select2.resid).keys())
        df_sasa = pd.DataFrame()

        signal.signal(signal.SIGINT, handler)
        logging.info(f'Running on pdb # {pdb}')
        select1 = atomsel(selection1)
        if (bp):
            sasa_l = []
            res_list = []
            #till = int(len(resids_select2) / 2) + 1
            for j, resids in enumerate(resids_select2):
                logging.info(f'{resids} in {selection1} and  {resids_select2[-j - 1]} in {selection2}')
                res_list.append(f'fc{resids}_sc{resids_select2[-j - 1]}')
                rest_resid = atomsel(f'({selection1} and resid \"{resids}\") or ({selection2} and resid \"{resids_select2[-j - 1]}\")')
                sasa = select1.sasa(srad=srad, restrict=rest_resid)
                sasa_l.append(sasa)
                signal.signal(signal.SIGINT, handler)
            df_sasa["resid_bp"] = res_list
            df_sasa[f'sasa_{time}'] = sasa_l
            logging.info(f'writing to file {sel1}_{sel2}_t{time}_bp.csv')
            df_sasa.to_csv(f'{pathx}/{sel1}_{sel2}_t{time}_bp.csv', index=None)

        else:
            sasa_l = []
            for resid in resids_select2:
                rest_resid =  atomsel(f'{selection2} and resid \"{resid}\"')
                sasa = select1.sasa(srad=srad,restrict=rest_resid)
                sasa_l.append(sasa)
                signal.signal(signal.SIGINT, handler)
            df_sasa["resid"] = resids_select2
            df_sasa[f'sasa_{time}'] = sasa_l
            logging.info(f'writing to file {sel1}_{sel2}_t{time}_notbp.csv')
            df_sasa.to_csv(f'{pathx}/{sel1}_{sel2}_t{time}_notbp.csv', index=None)



class comdist():
    def com_dist_tail(pdb,selection1,selection2,sec):
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        def conv(x):
            return str(f'\"{str(x)}\"')

        path = os.getcwd()
        os.chdir(path)
        if not os.path.exists("05_distance"):
            os.makedirs("05_distance")
        pathx = path + "/05_distance"
        logging.info(f'Starting..')

        mol=vmd.molecule.load('pdb',pdb)
        num_frame=vmd.molecule.numframes(mol)

        select1 = atomsel(f'{selection1}', frame=1)
        resids_select1 = list(OrderedDict.fromkeys(select1.resid).keys())

        select2 = atomsel(f'{selection2}', frame=1)
        resids_select2 = list(OrderedDict.fromkeys(select2.resid).keys())


        toappend=[[0,0,0]] # get first frame

        for i in range(num_frame):
            #signal.signal(signal.SIGINT, handler)
            logging.info(f'Running on frame # {i + 1}')

            #find a solution for minus signed nucleotide numbers
            strand2 = atomsel(f'({selection1} and resid {" ".join(map(conv,resids_select1[-6:]))}) or ({selection2} and '
                              f'resid {" ".join(map(conv,resids_select2[:6]))}) ',frame=i+1)
            strand1 = atomsel(f'({selection2} and resid {" ".join(map(conv,resids_select2[-6:]))}) or ({selection1} and '
                              f'resid {" ".join(map(conv,resids_select1[:6]))}) ',frame=i+1)

            strand1_c=strand1.center() #dna
            strand2_c=strand2.center() #dna


            histone = atomsel(f'({sec}) and helix and alpha',frame=i+1)
            histone_c = histone.center()

            strand1__histone=((histone_c[0] - strand1_c[0])**2 +  (histone_c[1] - strand1_c[1])**2 + \
                              (histone_c[2] - strand1_c[2])**2)**(1/2)
            strand2__histone=((histone_c[0] - strand2_c[0])**2 +  (histone_c[1] - strand2_c[1])**2 + \
                              (histone_c[2] - strand2_c[2])**2)**(1/2)
            toappend.append([i+1, strand1__histone, strand2__histone])  # append to others

        df_dist=pd.DataFrame(toappend,columns=["time","strand1","strand2"])
        df_dist = df_dist.iloc[1:]
        logging.info(f'Writing ..')
        df_dist.to_csv(pathx+"/distance.csv",index=None)

    def combine_com_dist_perres():

        pathx = os.getcwd()
        pathx_toreadz =  f'{pathx}/05_distance_perres/'

        if not os.path.exists(f'{pathx}/05_distance_perres_merged'):
            os.makedirs(f'{pathx}/05_distance_perres_merged', exist_ok=True)
        pathx_towrite = f'{pathx}/05_distance_perres_merged'

        df3 = pd.DataFrame()
        logging.info(f'Merging files ...')
        # get all files in a folder
        file_list = glob.glob(f'{pathx_toreadz}/*.csv')

        for file in file_list:
            df = pd.read_csv(f'{file}')
            df3 = df3.append(df)

        df3.to_csv(f'{pathx_towrite}/distance_perres.csv', index=None)


    def com_dist_perres(pdb,selection1,selection2,sec):
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        def conv(x):
            return str(f'\"{str(x)}\"')

        path = os.getcwd()
        os.chdir(path)
        if not os.path.exists("05_distance_perres"):
            os.makedirs("05_distance_perres", exist_ok=True)
        pathx = path + "/05_distance_perres"
        logging.info(f'Starting..')

        mol=vmd.molecule.load('pdb',pdb)
        time = pdb.split("md_")[1].split(".pdb")[0]
        sel1 = selection1.split(" ")[1]
        sel2 = selection2.split(" ")[1]

        select1 = atomsel(f'{selection1}')
        resids_select1 = list(OrderedDict.fromkeys(select1.resid).keys())

        select2 = atomsel(f'{selection2}')
        resids_select2 = list(OrderedDict.fromkeys(select2.resid).keys())

        logging.info(f'Running on #{time} pdb ..)')


        toappend=[[0,0,0,0,0]] # get first frame

        for i in range(len(resids_select2)):
            #signal.signal(signal.SIGINT, handler)

            strand2 = atomsel(f'({selection1} and resid {conv(resids_select1[-i-1])}) or ({selection2} and resid {conv(resids_select2[i])}) ')
            strand1 = atomsel(f'({selection2} and resid {conv(resids_select2[-i-1])}) or ({selection1} and resid {conv(resids_select1[i])})')


            strand1_c=strand1.center() #dna
            strand2_c=strand2.center() #dna


            histone = atomsel(f'({sec}) and helix and alpha')
            histone_c = histone.center()

            strand1__histone=((histone_c[0] - strand1_c[0])**2 +  (histone_c[1] - strand1_c[1])**2 + \
                              (histone_c[2] - strand1_c[2])**2)**(1/2)
            strand2__histone=((histone_c[0] - strand2_c[0])**2 +  (histone_c[1] - strand2_c[1])**2 + \
                              (histone_c[2] - strand2_c[2])**2)**(1/2)
            toappend.append([time,f'{resids_select1[-i-1]}_{resids_select2[i]}',f'{resids_select1[i]}_{resids_select2[-i-1]}', strand1__histone, strand2__histone])  # append to others

        df_dist=pd.DataFrame(toappend,columns=["time",f'strand1_chain{sel1}{sel2}',f'strand2_chain{sel2}{sel1}',"strand1_dist","strand2_dist"])
        df_dist = df_dist.iloc[1:]
        logging.info(f'Writing to {pathx} ..')
        df_dist.to_csv(f'{pathx}/distance_t{time}.csv',index=None)

    def com_dist_perres_merged(pdb,selection1,selection2,sec):
        def handler(sig, frame):
            raise KeyboardInterrupt("CTRL-C!")

        def conv(x):
            return str(f'\"{str(x)}\"')

        path = os.getcwd()
        os.chdir(path)
        if not os.path.exists("05_distance_perres"):
            os.makedirs("05_distance_perres", exist_ok=True)
        pathx = path + "/05_distance_perres"
        logging.info(f'Starting..')

        mol=vmd.molecule.load('pdb',pdb)
        time = pdb.split("md_")[1].split(".pdb")[0]
        sel1 = selection1.split(" ")[1]
        sel2 = selection2.split(" ")[1]

        select1 = atomsel(f'{selection1}')
        resids_select1 = list(OrderedDict.fromkeys(select1.resid).keys())

        select2 = atomsel(f'{selection2}')
        resids_select2 = list(OrderedDict.fromkeys(select2.resid).keys())

        logging.info(f'Running on #{time} pdb ..)')


        toappend=[[0,0,0,0,0]] # get first frame

        for i in range(len(resids_select2)):
            #signal.signal(signal.SIGINT, handler)

            strand2 = atomsel(f'({selection1} and resid {conv(resids_select1[-i-1])}) or ({selection2} and resid {conv(resids_select2[i])}) ',frame=i+1)
            strand1 = atomsel(f'({selection2} and resid {conv(resids_select2[-i-1])}) or ({selection1} and resid {conv(resids_select1[i])})',frame=i+1)


            strand1_c=strand1.center() #dna
            strand2_c=strand2.center() #dna


            histone = atomsel(f'({sec}) and helix and alpha',frame=i+1)
            histone_c = histone.center()

            strand1__histone=((histone_c[0] - strand1_c[0])**2 +  (histone_c[1] - strand1_c[1])**2 + \
                              (histone_c[2] - strand1_c[2])**2)**(1/2)
            strand2__histone=((histone_c[0] - strand2_c[0])**2 +  (histone_c[1] - strand2_c[1])**2 + \
                              (histone_c[2] - strand2_c[2])**2)**(1/2)
            toappend.append([time,f'{resids_select1[-i-1]}_{resids_select2[i]}',f'{resids_select1[i]}_{resids_select2[-i-1]}', strand1__histone, strand2__histone])  # append to others

        df_dist=pd.DataFrame(toappend,columns=["time",f'strand1_chain{sel1}{sel2}',f'strand2_chain{sel2}{sel1}',"strand1_dist","strand2_dist"])
        df_dist = df_dist.iloc[1:]
        logging.info(f'Writing to {pathx} ..')
        df_dist.to_csv(f'{pathx}/distance_t{time}.csv',index=None)




def split_pdbs(target):
    def handler(sig, frame):
        raise KeyboardInterrupt("CTRL-C!")

    path=os.getcwd()
    if not os.path.exists("02_pdbs"):
        os.makedirs("02_pdbs")
    pathxx = path + "/02_pdbs"

    tName = "md"
    fp = open(target, 'r')
    xx = []
    i = 0
    for line in fp:
        if line.startswith("ATOM"):
            xx.append(line.rstrip())
        elif line.startswith("TER"):
            xx.append(line.rstrip())
        elif line.startswith("HETATM"):
            xx.append(line.rstrip())
        elif line.startswith("END"):
            xx.append("END")
            aa = pd.DataFrame(xx)
            aa.to_csv(pathxx + "/" + tName + "_" + str(i) + ".pdb", index=None, header=False)
            xx = []
            i = i + 1
            logging.info(f'END stated {i} number of times')
        signal.signal(signal.SIGINT, handler)

    fp.close()





def clean():

    folders = list(os.walk("."))[1:]
    for dir, subdirs, files in folders:
        if len(os.listdir(dir)) == 0:
            logging.info(f'Deleting empty folder - {dir}')
            shutil.rmtree(dir)

    if os.path.exists('02_pdbs') and os.path.isdir('02_pdbs'):
        logging.info(f'Removing folder - 02_pdbs')
        shutil.rmtree('02_pdbs')

    if os.path.exists('09_nuc_e') and os.path.isdir('09_nuc_e'):
        logging.info(f'Removing folder - 09_nuc_e')
        shutil.rmtree('09_nuc_e')

    if os.path.exists('08_nuc_rc_sm') and os.path.isdir('08_nuc_rc_sm'):
        logging.info(f'Removing folder - 08_nuc_rc_sm')
        shutil.rmtree('08_nuc_rc_sm')

    if os.path.exists('07_nuc_com_xyz') and os.path.isdir('07_nuc_com_xyz'):
        logging.info(f'Removing folder - 07_nuc_com_xyz')
        shutil.rmtree('07_nuc_com_xyz')

    if os.path.exists('06_contact') and os.path.isdir('06_contact'):
        logging.info(f'Removing folder - 06_contact')
        shutil.rmtree('06_contact')


    if os.path.exists('04_sasa') and os.path.isdir('04_sasa'):
        logging.info(f'Removing folder - 04_sasa')
        shutil.rmtree('04_sasa')

    if os.path.exists('06_contact_perres') and os.path.isdir('04_sasa'):
        logging.info(f'Removing folder - 06_contact_perres')
        shutil.rmtree('06_contact_perres')













