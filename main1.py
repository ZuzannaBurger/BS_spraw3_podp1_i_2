from Bio.PDB.PDBList import PDBList
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
import os

pdbl = PDBList()
fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format="pdb")

pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure('4ywo', fetch_pdb)

ca_atoms = []
for model in structure:
    for chain in model:
        for res in chain:
            for atom in res.get_atoms():
                if atom.get_name() == 'CA':
                    ca_atoms.append(res['CA'])

    threshold = 8.0
    n = len(ca_atoms)
    contact_map = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            distance = np.linalg.norm(ca_atoms[i].coord - ca_atoms[j].coord)
            if distance <= threshold:
                contact_map[i, j] = 1
                contact_map[j, i] = 1

    plt.imshow(contact_map, cmap="binary", origin="upper")
    plt.title(f"Mapa kontaktów białka: {'4ywo'.upper()}")
    plt.xlabel("Indeks reszty")
    plt.ylabel("Indeks reszty")
    plt.colorbar(label="Kontakt (1=kontakt, 0=brak kontaktu)")
    plt.show()