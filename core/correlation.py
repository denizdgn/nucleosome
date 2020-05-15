
"""

Created on Sat April 25 2020
@author: Seyit Kale
@email: seyitkale@gmail.com

"""

import numpy as np
import os
from Bio.PDB.PDBParser import PDBParser
from sklearn.cluster import KMeans
# from sklearn.cluster import AgglomerativeClustering as agglo
from sklearn.metrics import normalized_mutual_info_score as mi
# from sklearn.metrics import mutual_info_score as mi
import logging
import sys

logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] %(message)s',
                    datefmt='%Y/%m/%d %H:%M:%S')




#
# Wrapper for BioPDB
#
def correlated(n_frames,n_bases,n_clusters,dna_chain_1,dna_chain_2,histone_chains):
    def get_residue_by_number(residues, number):
        """

        Handles Biophyton to regular array indexing

        """

        for residue in residues:
            if residue.id[1] == number:
                return residue


    #
    # Compute geometric center of a set of atoms
    #

    def get_base_center(residue):
        """

        Returns the geometric center of a set of atoms (a base)

        """

        residue_array = np.asarray(residue)

        com = np.zeros(3)
        n_atoms = len(residue_array)

        for i in range(0, n_atoms):
            res = residue_array[i].get_vector()
            for j in range(0, 3):
                com[j] += res[j]

        com /= n_atoms

        return com


    #
    # Compute the orientation vector of the base pair
    #


    def get_bp_vector(structure, index1, index2):
        """

        Returns the geometric center of a base pair

        """

        base1 = get_residue_by_number(structure[0][dna_chain_1], index1)
        base2 = get_residue_by_number(structure[0][dna_chain_2], index2)

        # The base vector yields higher MIs than the position (com) vector
        # norm vector performs worst.
        return get_base_center(base1) - get_base_center(base2)
        # return .5 * (get_base_center(base1) + get_base_center(base2))
        # return np.linalg.norm(get_base_center(base1) - get_base_center(base2))


    #
    # Return an array of base pair vectors in a structure
    #


    def get_bp_vectors(structure):
        """

        Return the geometric centers of all base pairs

        """

        nbp_1 = len(structure[0][dna_chain_1])
        com_1 = np.zeros([nbp_1, 3])

        counter = 0
        for residue in (structure[0][dna_chain_1]):
            com_1[counter] = get_base_center(residue)
            counter += 1

        nbp_2 = len(structure[0][dna_chain_2])
        com_2 = np.zeros([nbp_2, 3])

        counter = 0
        for residue in (structure[0][dna_chain_2]):
            com_2[counter] = get_base_center(residue)
            counter += 1

        if nbp_1 != nbp_2:
            print("WARNING: Mismatch in DNA strands")

        # print(str(nbp) + ' base pairs are detected.')

        base_pair_vectors = np.zeros([nbp_1, 3])

        for index in range(0, nbp_1):
            base_pair_vectors[index] = com_1[index] - com_2[nbp_1 - index - 1]

        return base_pair_vectors


    #
    # Parse one structure
    #


    # first_base = -72
    # last_base = 72
    # n_bases = last_base - first_base + 1

    parser = PDBParser(PERMISSIVE=1)


    def parse_structure(frame_number):
        """

        Parse a structure and return base-pair information

        """

        pathx = os.getcwd()
        folder_name = pathx+"/02_pdbs"

        pdb_filename = folder_name + "/md_" + str(frame_number) + ".pdb"

        return parser.get_structure("structure", pdb_filename)


    #
    # Declare arrays for residue information (COMs or vectors, all 3D)
    #

    bps_frames = np.zeros([n_bases, 3, n_frames])





    #
    # Go through each frame
    #

    for frame in range(0, n_frames):

        logging.info(f'Running on  {frame}..')

        structure = parse_structure(frame)

        bps = get_bp_vectors(structure)

        bps_frames[:, :, frame] = bps

    #
    # Go through all base pairs, cluster each to k states
    #

    #k_means = KMeans(n_clusters=n_clusters,n_jobs=-2)
    logging.info(f'Clustering..')
    k_means = KMeans(n_clusters=n_clusters)
    # The more sophisticated clustering algorithm does not perform better
    # than the simpler K-means.
    # k_means = agglo(n_clusters=n_clusters)

    dna_clusters = np.zeros([n_bases, n_frames])

    for base_pair in range(0, n_bases):
        dna_clusters[base_pair] = k_means.fit_predict(bps_frames[base_pair].T)

    # """
    mutual_information = np.zeros([n_bases, n_bases])

    for ibp1 in range(0, n_bases):
        for ibp2 in range(0, n_bases):
            bp1 = dna_clusters[ibp1]
            bp2 = dna_clusters[ibp2]

            mutual_information[ibp1][ibp2] = mi(bp1, bp2)

    pathx = os.getcwd()
    if not os.path.exists(pathx + "/10_mutual_information"):
        os.makedirs(pathx + "/10_mutual_information",exist_ok=True)
    pathxz = pathx + "/10_mutual_information"

    np.save(f'{pathxz}/mi', mutual_information)
    #np.savetxt(f'{pathxz}/mi.csv', mutual_information, delimiter=",")

    # """
