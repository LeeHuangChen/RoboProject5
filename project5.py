from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain
from MMTK.Collections import Collection
from MMTK.ForceFields import Amber94ForceField
#from MMTK.Visualization import view
from Scientific.Visualization import VMD; module = VMD

# Load the PDB file.
#configuration = PDBConfiguration('1nxb.pdb')
configuration = PDBConfiguration('110d.pdb')

# Construct the peptide chain objects.
#chains = configuration.createPeptideChains()

# Make the protein object.
#protein = Protein(chains)

# The protein contains only one chain
chain = configuration.peptide_chains[0]

# Number of residues in Protein
numberOfResidues = len(chain)

# Access all residues
#for residue in chain: do sth

# Keep the substructure formed by the residues 23–42
# p[i:j] yields the subchain of chain p from residue number i up to but excluding residue number j
sc = chain[23:43]

# Create a universe and put the chain in it
#universe = InfiniteUniverse()
universe = InfiniteUniverse(Amber99ForceField())
universe.addObject(sc)
# Current configuration of universe
universe.configuration()
# List of atoms in universe
universe.atomList()

#view(universe)

# Visualization
graphics = protein.graphicsObjects(graphics_module = module,
                                   model = 'backbone', color = 'red')
scene = VMD.Scene(graphics)
scene.view()

