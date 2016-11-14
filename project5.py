from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain
from MMTK.Proteins import Protein
from MMTK.Collections import Collection
from MMTK.ForceFields import Amber94ForceField
#from MMTK.Visualization import view
from Scientific.Visualization import VMD; module = VMD

# Load the PDB file.
#configuration = PDBConfiguration('1nxb.pdb')
#configuration = PDBConfiguration('110d.pdb')
#configuration = PDBConfiguration('insulin.pdb')
#configuration = PDBConfiguration('bALA1.pdb')
configuration = PDBConfiguration('2YCC_mod2.pdb')

# Construct the peptide chain objects.
chains = configuration.createPeptideChains()

# Make the protein object.
#protein = Protein(chains)

# The protein contains only one chain
#chain = configuration.peptide_chains[0]
chain = chains[0]
#sc = chain[11:27]

# Number of residues in Protein
numberOfResidues = len(chain)

# Access all residues
goalAngles = []
for i in range(1,len(chain)-1): 
	print i,	
	goalAngles.append(chain[i].chiAngle().getValue())
	chain[i].phiAngle().setValue(.5)
	chain[i].psiAngle().setValue(.5)


print goalAngles


# p[i:j] yields the subchain of chain p from residue number i up to but excluding residue number j
#sc = chain[23:43]

# access the i-th residue
#r1 = chain[1]

# DihedralAngle 
#angle = r1.chiAngle()

# Create a universe and put the chain in it
universe = InfiniteUniverse(Amber94ForceField())
protein = Protein(chain)
universe.addObject(protein)
energy = universe.energy()
print energy
# Current configuration of universe
universe.configuration()
# List of atoms in universe
universe.atomList()

#view(universe)
sc = chain[1:3]
# Visualization
graphics = sc.graphicsObjects(graphics_module = module,
                                   model = 'backbone', color = 'red')
scene = VMD.Scene(graphics)
scene.view()

