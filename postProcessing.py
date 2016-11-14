
try:
    from ompl import base as ob
    from ompl import geometric as og
    from ompl import util as ou
except:
    # if the ompl module is not in the PYTHONPATH assume it is installed in a
    # subdirectory of the parent directory called "py-bindings."
    from os.path import abspath, dirname, join
    import sys
    sys.path.insert(0, join(dirname(dirname(abspath(__file__))),'py-bindings'))
    
    from ompl import base as ob
    from ompl import geometric as og

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain
from MMTK.Proteins import Protein
from MMTK.Collections import Collection
from MMTK.ForceFields import Amber94ForceField

numberOfResidues = 16
n=numberOfResidues*2
space = ob.CompoundStateSpace()

for i in range(n): space.addSubspace(ob.SO2StateSpace(),1)

ss = og.SimpleSetup(space)
   
pds=ob.PlannerDataStorage()
data=ob.PlannerData(ss.getSpaceInformation())
pds=ob.PlannerDataStorage()
pds.load("PlannerData.txt", data)
#print data.numStartVertices()
states = data.extractStateStorage()
#data.toBoostGraph()
print states[1].value





