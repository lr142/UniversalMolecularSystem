import sys
sys.path.append("..")
from Utility import *
from UniversalMolecularSystem import *
from MOL2File import *


def BreakupMoleculeByConnectivity(molecule):
    # returns a MolecularStructure
    # Make sure that it has the bonding map 'bondedTo' properly constructed.

    atomCount = len(molecule.atoms)

    # Initially, every atom is associated with a parent (itself) and an empty set of children
    parent = [i for i in range(atomCount)]
    children = [set() for i in range(atomCount)]

    # Merge two 'flat' trees
    def TreeMerge(n1,n2):
        p1 = parent[n1]
        p2 = parent[n2]
        if p1 == p2:
            return
        for c in children[p2]:
            parent[c] = p1
            children[p1].add(c)

        children[p1].add(p2)
        parent[p2] = p1
        children[p2] = set()


    # Now we check the bonding map
    for index, connected in enumerate(molecule.bondedTo):
        # For each atom
        for toAtomIndex in connected:
            # Only need to check its connectivity with previous atoms
            if toAtomIndex < index:
                TreeMerge(toAtomIndex,index)

    # Debugging, Display the tree
    # for i in range(atomCount):
    #     if parent[i] == i:  # This is a root node
    #         output("Root : {}, with {} children: \n{}".format(i,len(children[i]),children[i]))
    



def TestBreakupMolecule():
    m = MolecularSystem()
    m.Read(MOL2File(),'Coal70A.mol2')
    import random
    random.shuffle(m.molecules[0].atoms)
    BreakupMoleculeByConnectivity(m.molecules[0])
    from XYZFile import XYZFile

    # file = open('dump.xyz','w')
    # output.setoutput(file)
    # m.Write(XYZFile())
    # output.setoutput(None)
    # file.close()

TestBreakupMolecule()


