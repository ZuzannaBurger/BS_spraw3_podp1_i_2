from Bio.PDB import PDBList, PDBParser, vectors
import seaborn as sns
import matplotlib.pyplot as plt

pdbl = PDBList()
fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')


parser = PDBParser(QUIET=True)
structure = parser.get_structure('4YWO', fetch_pdb)

def calculate_phi_psi(structure):
    phi_angles = []
    psi_angles = []

    for model in structure:
        for chain in model:
            residues = list(chain)
            for i in range(1, len(residues) - 1):
                try:

                    c_prev = residues[i - 1]['C'].get_vector()
                    n_curr = residues[i]['N'].get_vector()
                    ca_curr = residues[i]['CA'].get_vector()
                    c_curr = residues[i]['C'].get_vector()
                    n_next = residues[i + 1]['N'].get_vector()

                    phi = vectors.calc_dihedral(c_prev, n_curr, ca_curr, c_curr)
                    phi_angles.append(phi * (180.0 / 3.141592653589793))

                    psi = vectors.calc_dihedral(n_curr, ca_curr, c_curr, n_next)
                    psi_angles.append(psi * (180.0 / 3.141592653589793))

                except KeyError:
                    continue

    return phi_angles, psi_angles


phi_angles, psi_angles = calculate_phi_psi(structure)

phi_angles = [phi for phi in phi_angles if -180 <= phi <= 180]
psi_angles = [psi for psi in psi_angles if -180 <= psi <= 180]

plt.figure(figsize=(10, 8))
sns.kdeplot(x=phi_angles, y=psi_angles, cmap="BuPu", fill=True, levels=50, thresh=0.05)
sns.scatterplot(x=phi_angles, y=psi_angles, color='black', s=5, alpha=0.5)

plt.title('Ramachandran Diagram for Protein 4YWO', fontsize=16)
plt.xlabel('Phi (°)', fontsize=14)
plt.ylabel('Psi (°)', fontsize=14)
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.grid(True, linestyle='--', alpha=0.5)
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.show()
