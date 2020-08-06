import sys
sys.path.append("..")
from Utility import *
from UniversalMolecularSystem import *
from MOL2File import *
from XYZFile import *


def BreakupMolecule(originalMolecule, atomBelongToNewMoleculeSerialMap):
    # This function breaks up a molecule and return a MolecularSystem.
    # If successful, a NEWLY CREATED MolecularSystem, along with newly created constituent Molecules
    #    and Atoms shall be returned. Otherwise it returns None.
    # The breakup scheme is 'atomBelongToNewMoleculeSerialMap', which maps each atom serial to a new molecule serial
    # Requirements:
    # 1. atom serials in the original molecule must be unique
    # 2. each atom must belong to one and only one molecule in the new system
    # Features:
    # 1. After the breakup, some bonds may become intermolecular, which are move in the 'interMolecularBonds' property
    #    in the new MolecularSystem. Those bond that are entirely within a single molecule in the new system becomes a
    #    part of the newly created molecule.
    # 2. While the serials of new molecules are specified by the user, all new molecules will duplicate the name of the original one
    # 3. Molecules in the new system are sorted according to the smallest index of its constituent atoms in the original molecule.
    #    For example, if atom [2 3 6],[1 4], and [5] are assigned to molecule x, y, and z in the new system, then molecule y (atom [1 4])
    #    shall come first, then molecule x ([2,3,6]) and molecule z ([5]).


    foundSerial = set()
    newMoleculeSerialToIndexMap = {}
    for atom in originalMolecule.atoms:
        # Check that atom serials are unique, and each atom is assigned to a new molecule
        if atom.serial in foundSerial:
            error("BreakupMolecule() fails because atom serials are not unique. Serial {} appears at least twice.".format(atom.serial),False)
            return None
        else:
            foundSerial.add(atom.serial)

        if not atom.serial in atomBelongToNewMoleculeSerialMap:
            error("BreakupMolecule() fails because the atom with serial {} are not assigned to a new molecule.".format(atom.serial),False)
            return None
        else:
            molSerial = atomBelongToNewMoleculeSerialMap[atom.serial]
            if not molSerial in newMoleculeSerialToIndexMap:
                molIndex = len(newMoleculeSerialToIndexMap)
                newMoleculeSerialToIndexMap[molSerial] = molIndex

    molCount = len(newMoleculeSerialToIndexMap)   # number of new molecules found in the first scan
    if molCount == 0 or len(originalMolecule.atoms) == 0:  # Do nothing upon empty molecule or empty splitting scheme
        return None

    newMS = MolecularSystem()                     # creates new Molecule objects and set their properties
    newMS.molecules = [Molecule() for i in range(molCount)]
    for molSerial in newMoleculeSerialToIndexMap:
        index = newMoleculeSerialToIndexMap[molSerial]
        mol = newMS.molecules[index]
        mol.serial = molSerial
        mol.name = originalMolecule.name
        mol.type = originalMolecule.type

    def FindMoleculeFromAtomSerial(atomSerial):
        molSerial = atomBelongToNewMoleculeSerialMap[atomSerial]
        molIndex = newMoleculeSerialToIndexMap[molSerial]
        return newMS.molecules[molIndex]

    import copy
    # copy atoms into the new system
    serialGlobalToLocalMap = {}
    for atom in originalMolecule.atoms:
        newMol = FindMoleculeFromAtomSerial(atom.serial)
        newAtom = copy.copy(atom)               # A 'shallow copy' shall do
        newAtom.systemwideSerial = newAtom.serial   # Sets the systemwide serial.
        # Also need to renumber atoms within each molecule
        newAtom.serial = "{}".format(len(newMol.atoms)+1)
        serialGlobalToLocalMap[newAtom.systemwideSerial] = newAtom.serial

        newMol.atoms.append(newAtom)

    # copy bonds into the new system. Need to distinguish between inter- and intra-molecular bonds
    for bond in originalMolecule.bonds:
        mol1 = FindMoleculeFromAtomSerial(bond.atom1)
        mol2 = FindMoleculeFromAtomSerial(bond.atom2)
        if mol1 == mol2:    # comparing two references, shall be fine  # intra-molecular
            # In this case, must switch to the local (within a molecule) serial
            newBond = copy.copy(bond)
            newBond.atom1 = serialGlobalToLocalMap[bond.atom1]
            newBond.atom2 = serialGlobalToLocalMap[bond.atom2]
            mol1.bonds.append(newBond)
        else:  # intermolecular
            if newMS.interMolecularBonds == None:
                newMS.interMolecularBonds = []
            newMS.interMolecularBonds.append(copy.copy(bond))

    # for m in newMS.molecules:
    #     error.turn_on()
    #     m.CheckConsistency()
    # exit()
    return newMS

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
    bondedTo = molecule.BondedMap()
    for index, connected in enumerate(bondedTo):
        # For each atom
        for toAtomIndex in connected:
            # Only need to check its connectivity with previous atoms
            if toAtomIndex < index:
                TreeMerge(toAtomIndex,index)

    #Debugging, Display the tree
    # ms = MolecularSystem()
    # for i in range(atomCount):
    #     if parent[i] == i:  # This is a root node
    #         ms.molecules.append(Molecule())
    #         output("Root : {}, with {} children: \n{}".format(i,len(children[i]),children[i]))
    #         ms.molecules[-1].atoms.append(molecule.atoms[i])
    #         for c in children[i]:
    #             ms.molecules[-1].atoms.append(molecule.atoms[c])
    # output.setoutput(open('dump.xyz','w'))
    # ms.Write(XYZFile())
    # exit()

    # Create the belongTo map
    belongTo = {}
    moleculeSerialCounter = 0
    # 1st round, scan for all root nodes and create molecules as the scan goes
    for i, atom in enumerate(molecule.atoms):
        if parent[i] == i: # Root node found, which suggests a new molecule should be created
            moleculeSerialCounter += 1
            belongTo[atom.serial] = '{}'.format(moleculeSerialCounter)
    # 2nd round, sets the belongTo map properly.
    for i, atom in enumerate(molecule.atoms):
        # Set its belongTo molecule to its parent's. For root nodes, this operation has no effect
        parentAtom = molecule.atoms[parent[i]]
        belongTo[atom.serial] = belongTo[parentAtom.serial]

    # All set, break it up!
    return BreakupMolecule(molecule,belongTo)

def BreakupMoleculeRandomly(molecule,numberOfMolecules):
    belongTo = {}
    import random
    for a in molecule.atoms:
        belongTo[a.serial] = random.randint(1,numberOfMolecules)
    return BreakupMolecule(molecule,belongTo)

def ReduceSystemToOneMolecule(originalSystem):
    # This function combines all molecules in a system into a single molecule.
    # If successful, the NEWLY CREATED molecule shall be returned. Otherwise it returns None
    # Requirement:
    # 1. If there are intermolecular bonds, each atom must have a unique system-wide serial number
    # Features:
    # 1. If each atom has a systemwide serial, this serial will replace its original, molecular-wide serial.
    #    If no systemwide serial is available, a consecutive number will be assigned to each atom.
    # 2. Name and type of the combined molecule is set to be the same as the first molecule.

    if len(originalSystem.molecules) == 0:    # Quit on emtpy systems:
        return None

    systemwideSerialsSet = set()
    systemwideSerialsOK = False
    if originalSystem.interMolecularBonds != None and len(originalSystem.interMolecularBonds) > 0:
        # Check that each atom has a unique systemwide serial
        for m in originalSystem.molecules:
            for a in m.atoms:
                if a.systemwideSerial == None or a.systemwideSerial in systemwideSerialsSet:
                    error("ReduceSystemToOneMolecule() fails because the some atom doesn't have unique systemwideSerial"
                          " when intermolecular bonds are present",False)
                    return None
                else:
                    systemwideSerialsSet.add(a.systemwideSerial)
        systemwideSerialsOK = True


    newMolecule = Molecule()
    newMolecule.name = originalSystem.molecules[0].name
    newMolecule.type = originalSystem.molecules[0].type
    import copy
    for m in originalSystem.molecules:

        localSerialToGlobalSerialMap = {}
        for a in m.atoms:
            newAtom = copy.copy(a)
            if not systemwideSerialsOK: # If there is originally no systemwide serial, assign one to it
                newAtom.systemwideSerial = "{}".format(len(newMolecule.atoms)+1)
            localSerialToGlobalSerialMap[newAtom.serial] = newAtom.systemwideSerial
            newAtom.serial = newAtom.systemwideSerial
            newMolecule.atoms.append(newAtom)
        for b in m.bonds:
            newBond = copy.copy(b)
            try:
                newBond.atom1 = localSerialToGlobalSerialMap[newBond.atom1]
                newBond.atom2 = localSerialToGlobalSerialMap[newBond.atom2]
            except:
                error("Can't find atom {} or atom {}!".format(newBond.atom1,newBond.atom2))
            newMolecule.bonds.append(newBond)

    if originalSystem.interMolecularBonds != None:
        for b in originalSystem.interMolecularBonds:
            newBond = copy.copy(b)
            newMolecule.bonds.append(newBond)

    return newMolecule

def TestBreakupMolecule():
    ms = MolecularSystem()
    ms.Read(MOL2File(),'Coal70A.mol2')
    ms.molecules[0].Summary()

    ms_split = None;
    for i in range(2):
        output("Round {}".format(i))
        import random
        random.shuffle(ms.molecules[0].atoms)
        ms_split = BreakupMoleculeRandomly(ms.molecules[0],random.randint(1,200))
        #ms_split = BreakupMoleculeByConnectivity(ms.molecules[0])
        ms_split.Summary()
        newMol = ReduceSystemToOneMolecule(ms_split)

        newMol.Summary()
        ms = MolecularSystem()
        ms.molecules = [newMol]

    output.setoutput(open('dump.mol2','w'))
    ms_split = BreakupMoleculeByConnectivity(ms.molecules[0])
    ms_split.Write(MOL2File())



TestBreakupMolecule()


