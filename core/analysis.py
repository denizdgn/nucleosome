import os
import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool
import vmd
import sys
from vmd import *
from collections import OrderedDict
import logging
import signal
import pathos
import glob
pool= Pool(pathos.multiprocessing.cpu_count()-2)


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')


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
        for file in files_csv:
            df = pd.read_csv(file)
            df3[df.columns[1]] = df.iloc[:,1]
            signal.signal(signal.SIGINT, handler)
        df3["mean"] = df3.iloc[:, 1:].mean(axis=1)
        df4=pd.concat([pd.DataFrame(df.iloc[:,0]),df3],axis=1)
        df4.to_csv(f'{pathx_towrite}/{folder}.csv', index=None)

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
        df_res_cont_num.to_csv(pathx + "/contact_from_" + sel1 + "_to_" + sel2 + f'_t{time}.csv',
                                   index=None)













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




#def sasa(pdb:str,selection1:str,restrict:str,srad:float):
#
#    path=os.getcwd()
#    os.chdir(path)
#    if not os.path.exists("04_sasa"):
#        os.makedirs("04_sasa")
#    pathx = path + "/04_sasa"
#    mol = vmd.molecule.load('pdb', pdb)
#
#    num_frame = vmd.molecule.numframes(mol)
#    df_sasa = pd.DataFrame(columns=["sasa", "time"])
#    sasa_l = []
#    for i in range(num_frame):
#        logging.info(f'Running on frame # {i + 1}')
#        select1 = atomsel(selection1, frame=i + 1)
#        rest1 = atomsel(restrict, frame=i + 1)
#        sasa = select1.sasa(srad=srad, restrict=rest1)
#        sasa_l.append(sasa)
#
#    df_sasa["time"] = range(1,num_frame+1)
#    df_sasa["sasa"] = sasa_l
#    sel1 = selection1.replace(" ", "")
#    df_sasa.to_csv(pathx + "/sasa_" + sel1 +"_"+restrict +".csv", index=None)




def sasa(pdb:str,selection1:str,selection2:str,srad:float,bp=False):
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
    os.chdir(path)

    if not os.path.exists("04_sasa"):
        os.makedirs("04_sasa")
    pathx = path + "/04_sasa"
    logging.info(f'Reading trajectory .. ')
    mol = vmd.molecule.load('pdb', pdb)
    num_frame=vmd.molecule.numframes(mol)

    select1 = atomsel(f'{selection1}', frame=1)
    select2 = atomsel(f'{selection2}', frame=1)

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

    #num_frame = 2
    df_sasa = pd.DataFrame()
    for i in range(num_frame):
        signal.signal(signal.SIGINT, handler)

        logging.info(f'Running on frame # {i + 1}')
        select1 = atomsel(selection1, frame=i + 1)

        if (bp):
            sasa_l = []
            res_list = []
            #till = int(len(resids_select2) / 2) + 1
            for j, resids in enumerate(resids_select2):
                logging.info(f'{resids} in {selection1} and  {resids_select2[-j - 1]} in {selection2}')
                res_list.append(f'fc{resids}_sc{resids_select2[-j - 1]}')
                rest_resid = atomsel(f'({selection1} and resid \"{resids}\") or ({selection2} and resid \"{resids_select2[-j - 1]}\")', frame=i + 1)
                sasa = select1.sasa(srad=srad, restrict=rest_resid)
                sasa_l.append(sasa)
            df_sasa["resid_bp"] = res_list
            df_sasa[f'sasa_{i + 1}'] = sasa_l

        else:
            sasa_l = []
            for resid in resids_select2:
                #print(resid)
                rest_resid =  atomsel(f'{selection2} and resid \"{resid}\"', frame=i + 1)
                sasa = select1.sasa(srad=srad,restrict=rest_resid)
                sasa_l.append(sasa)
            df_sasa["resid"] = resids_select2
            df_sasa[f'sasa_{i+1}'] = sasa_l


    df_sasa["mean"]=df_sasa.iloc[:,1:].mean(axis=1)
    sel2=selection2.replace(" ", "")

    if (bp):
        df_sasa.to_csv(pathx+"/sasa_"+sel2+"_bp.csv",index=None)
    elif (not bp):
        df_sasa.to_csv(pathx + "/sasa_" + sel2 + ".csv", index=None)





def com_dist(pdb,selection1,selection2,sec):
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
    #def multiple(i,selection1,selection2,)
    for i in range(num_frame):
        signal.signal(signal.SIGINT, handler)
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



















