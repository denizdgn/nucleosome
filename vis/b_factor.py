import pandas as pd
import io
from collections import OrderedDict

def read_pdb(target):
    """
    :param target: single pdb file
    :return:
    """

    # path = os.getcwd()
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
    df_filt.drop(["temp_fac", "segid", "element_sym", "charge"], axis=1, inplace=True)
    df_filt = df_filt.fillna('')

    return df_filt


def generate_bfactor_per_atom(trajectory, bfactor_per_residue):
    """

    :param trajectory:
    :param bfactor_per_residue:
    :return:
    """

    bfactor_per_atom = np.zeros(trajectory.shape[0])
    residues = trajectory.apply(lambda x: x["chainID"] + "_" + x["res_name"] + "_"+ str(int(x["res_num"])), axis=1)


    for residue in range(len(residues)):
        chain = pdb[pdb.chainID == residues[0].split("_")[0]]
        residue_atoms = chain[chain.res_num == int(residues[0].split("_")[2])]

        residue_bfactor = bfactor_per_residue[residue]

        for atom in range(len(residue_atoms)):
            residue_atoms.occupancy[atom] = residue_bfactor

    return bfactor_per_atom




pathx="/Volumes/DENIZ/kalelab/Project_biowulf/MDx/round_1/cenpa/analysisNuc/"
os.chdir(f'{pathx}')
os.makedirs(f'{pathx}/00_graph/bfactor/', exist_ok=True)

pdb = read_pdb(f'{pathx}/02_pdbs/md_0.pdb')
pdb.res_num = pdb.res_num.astype(int)


first_bp=mi[0,]


b_factor = generate_bfactor_per_atom(pdb, first_bp)





generate_bfactor_per_atom(pdb,first_bp)
