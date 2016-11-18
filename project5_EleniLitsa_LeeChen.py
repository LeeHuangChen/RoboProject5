#!/usr/bin/env python


# Author: Lee Chen, Eleni Litsa

#import the modules
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
#from MMTK.Visualization import view
from Scientific.Visualization import VMD; module = VMD

import numpy as np

#global variables

#protein subchains
sc=None
sc2=None


numberOfResidues=None  #number of residues in the subchain

#containers for the proteins
universe=None
universe2=None

#two variables for the proteins
chain=None
chain2=None

#Threshold energy for the state validator
eThresh = None

#amount of time to conduct the calculation
time=None

# Custom projection class for KPIECE1
class MyProjectionEvaluator(ob.ProjectionEvaluator):
    def __init__(self, space, cellSizes):
        super(MyProjectionEvaluator, self).__init__(space)
        #self.setCellSizes(cellSizes)

    def getDimension(self):
        return 2

    def project(self, state, projection):
        projection[0] = state[0].value
        projection[1] = state[1].value
# Custom Motion Validator 
# Ran into python compilation problem while using this class
# This program uses instead a state validator for the planning
class MyMotionValidator(ob.MotionValidator):
    def __init__(self,spaceinfo):
        super(MyMotionValidator, self).__init__(spaceinfo)
    def checkMotion(s1, s2):
        avogadro=6.0221409e23
        kb = 1.3806488e-23
        kb_corr=kb/1000
        T=300.0
        
        j=0
        for i in range(0,len(sc)):
            sc[i].phiAngle().setValue(s1[j].value)
            j=j+1

            sc[i].psiAngle().setValue(s1[j].value)
            j=j+1
        einit=universe.energy()

        j=0
        for i in range(0,len(sc)):
            sc[i].phiAngle().setValue(s2[j].value)
            j=j+1

            sc[i].psiAngle().setValue(s2[j].value)
            j=j+1
        enext=universe.energy()

        exponent=(1.0*einit-enext)/(kb_corr*T*avogadro)
        #print "(",einit,",",enext,")",

        threshold=np.exp(exponent)

        rng=ou.RNG()
        if(enext<=einit):
            #print ""
            return True
        else:
            #print ":(",threshold,",",exponent,")"
            rValue=rng.uniform01()
            if(rValue<threshold):
                einit=enext
                return True
            else:
                return False
    def checkMotion(s1, s2, lastValid):
        lastValid=[]
        checkMotion(s1,s2)

# A state validator that uses the last accepted high energy state
# as the threshold to accept or reject a new state.
def isStateValid(state):
    global eThresh
    j=0
    avogadro=6.0221409e23
    kb = 1.3806488e-23
    kb_corr=kb/1000
    
    #print kb
    T=300.0
    #T=300e22
    for i in range(0,len(sc)):
        sc2[i].phiAngle().setValue(state[j].value)
        j=j+1

        sc2[i].psiAngle().setValue(state[j].value)
        j=j+1
    enext=universe2.energy()
    
    exponent=(1.0*eThresh-enext)/(kb_corr*T*avogadro)

    
    threshold=np.exp(exponent)
    
    #print "(",eThresh,",",enext,")",

    rng=ou.RNG()
    if(enext<=eThresh):
        #print ""
        return True
    else:
        #print ":(",threshold,",",exponent,")"
        rValue=rng.uniform01()
        if(rValue<threshold):
            eThresh=enext
            return True
        else:
            return False


def planWithSimpleSetup():
    global eThresh
    n=numberOfResidues*2

    #Construct the state space
    space = ob.CompoundStateSpace()
    for i in range(0,n):
        space.addSubspace(ob.SO2StateSpace(),1)

    #setup projection for the KPIECE planner
    myProjection = MyProjectionEvaluator(space,0.1)
    space.registerDefaultProjection(myProjection)
    


    # create a simple setup object
    ss = og.SimpleSetup(space)
    #add a state validity checker
    ss.setStateValidityChecker(ob.StateValidityCheckerFn(isStateValid))
    #ran into python problem while trying to use a custom motion validator
    #ss.getSpaceInformation().setMotionValidator(MyMotionValidator(ss.getSpaceInformation()))
    

    start = ob.State(space)
    goal = ob.State(space)

    # Access all residues and set the start and goal states
    # The goal state is defined with 0 for all angles, which
    # should be a very high energy state.
    j=0
    for i in range(0,len(sc)): 
        start[j]=sc[i].phiAngle().getValue()
        goal[j]=0
        j=j+1

        start[j]=sc[i].psiAngle().getValue()
        goal[j]=0
        j=j+1


    # Add both proteins into their on "containers" (i.e. universes)
    global universe 
    universe = InfiniteUniverse(Amber94ForceField())
    protein = Protein(chain)
    universe.addObject(protein)
    energy = universe.energy()
    global universe2
    universe2 = InfiniteUniverse(Amber94ForceField())
    protein2=Protein(chain2)
    universe2.addObject(protein2)

    #prints the initial energy for testing
    #print energy

    #set the initial energy as the current threshold energy
    eThresh=universe.energy()
    ss.setStartAndGoalStates(start, goal)

    
    #set planner
    planner=og.KPIECE1(ss.getSpaceInformation())
    #planner=og.RRT(ss.getSpaceInformation())
    ss.setPlanner(planner)
    ss.setup()
    planner.setGoalBias(0)

    #set the range of the planner to a very small value
    #this is because the protein has many high energy conformations
    prange=planner.getRange()
    const=1000
    prangeNew=prange/const
    planner.setRange(prangeNew)
    #set the goal bias to 0 so we can explore the configuration space
    #without goal biasing
    planner.setGoalBias(0)
    #prints the new and old range for testing
    #print "range:(",prange,",",prangeNew,")"

    global time
    solved = ss.solve(time)

    if solved:
        # try to shorten the path
        ss.simplifySolution()
        # print the simplified path
        print(ss.getSolutionPath())
        
    pds=ob.PlannerDataStorage()
    plannerdata=ob.PlannerData(ss.getSpaceInformation())
    planner.getPlannerData(plannerdata)
    pds.store(plannerdata,"PlannerData.txt")
    


if __name__ == "__main__":
    #load the first copy of the protein we are studying
    configuration = PDBConfiguration('2YCC_mod3.pdb')
    chains = configuration.createPeptideChains()
    chain = chains[0]
    sc = chain[11:20]

    #load the second copy of the protein we are studying
    configuration2 = PDBConfiguration('2YCC_mod3.pdb')
    chains2 = configuration.createPeptideChains()
    chain2 = chains2[0]
    sc2 = chain2[11:20]

    #ask the user on how long the experiment should run
    seconds = float(raw_input("How many seconds do you want to run the calculation? "))
    minutes = float(raw_input("How many minutes do you want to run the calculation? "))
    hours = float(raw_input("How many hours do you want to run the calculation? "))
    time=seconds+60*(minutes)+60*60*hours

    numberOfResidues = len(sc)
    planWithSimpleSetup()
