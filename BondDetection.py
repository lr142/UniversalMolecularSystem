# Automatically detect chemical bonds based on element types and distances.
# It should be implemented in a more advanced way, for example through machine learning...
# I'll work on this later...

# The class BondRules and the function MainAsXYZ2Mol2 are from legacy code, just for a reference.
# The class Bond is an older version of the Bond Class, which should also be rewritten.

# The [bond_rule_file] is a txt files (in csv format) which defines bond rules, they can include multiple lines,
# each line is like below:
# AtomType1 AtomType2 MinBondLength MaxBondLength BondType
# For example in the following lines
# C,O,1.2,1.6,1
# C,O,1.1,1.2,ar
# The above two lines ask the program to judge any two C, O atoms that are within 1.2 Angstrom(A) but at least 1.1 A
# apart form a aromatic bond, any two C, O atoms that are at least 1.2 A apart but within 1.6 A from each other form
# a single bond. Rules are checked in reverse order, which means that a rule appears later in the file overrides rules
# that appeared earlier, if it contradicts with the earlier rules.
# Other specifications include:
# 1. In each rule, the order of names of two atoms doesn't matter, meaning "C O 1.2 1.6 1" and "O C 1.2 1.6 1" are
#    exactly the same.
# 2. The 3rd parameter of this program, [bond_rule_file], is optional. There is a default bond rule file, named
#    "xyz2mol2.default.bond.rules", comes together with the xyz2mol2.py program and read by the program each time
#    it runs, to recognize some commonly encountered bonds in organic molecules. However a user defined [bond_rule_file]
#    can be supplied if some bonds are not recognized or special rules apply. Remember that rules in [bond_rule_file]
#    will override the default rules, if there are contradictions. Futhermore, rules appear later supersede previous one.
# 3. The 5th segment of each line in bond_rule_file defines the bond type, it can be any string that conforms with the
#    TRIPOS mol2 file format specification, such as 1, 2, 3, ar, am, etc. More information can be obtained from
#    www.tripos.com.

import sys
from Utility import *
from UniversalMolecularSystem import *
from XYZFile import *
from MOL2File import *
from MolecularManipulation import *
from PDBFile import *

class NeighborList:
    minx:float
    miny:float
    minz:float
    maxx:float
    maxy:float
    maxz:float
    gridSize:float
    Nx: int
    Ny: int
    Nz: int
    neighborList: [set]
    atomList: [Atom]
    #grid: a 3-dimensional array of sets, recording which atoms are in each grid
    def __init__(self, listOfAtoms, gridSize):

        if len(listOfAtoms) == 0:
            error("In NeighborList.__init__(), listOfAtoms must not be empty")
        if gridSize < 0:
            error("In NeighborList.__init__(), gridSize = {} must be positive".format(gridSize))

        self.minx = self.miny = self.minz = 1E10
        self.maxx = self.maxy = self.maxz = -1E10
        self.gridSize = float(gridSize)
        self.atomList = listOfAtoms
        for a in listOfAtoms:
            self.minx = min(a.x,self.minx)
            self.miny = min(a.y,self.miny)
            self.minz = min(a.z,self.minz)
            self.maxx = max(a.x,self.maxx)
            self.maxy = max(a.y,self.maxy)
            self.maxz = max(a.z,self.maxz)
        from math import floor
        self.Nx = int(floor( (self.maxx - self.minx) / self.gridSize)) + 1
        self.Ny = int(floor( (self.maxy - self.miny) / self.gridSize)) + 1
        self.Nz = int(floor( (self.maxz - self.minz) / self.gridSize)) + 1
        self.grid = [[[ set() for k in range(self.Nz)] for j in range(self.Ny) ] for i in range(self.Nx) ]

        # 2nd pass, assign atoms to each grid:
        for a in listOfAtoms:
            ix,iy,iz = self._find_grid_(a.x,a.y,a.z)
            self.grid[ix][iy][iz].add(a)

        # 3rd pass, create the neighborList
        N = len(listOfAtoms)
        self.neighborList = [set() for i in range(N)]

        def allAtomsInNeighborsOfAGridBlock(ix,iy,iz):
            # Find all atoms in the 9 blocks (including itself) around grid[ix, iy, iz]
            result = set()
            for i in range(max(ix-1,0),min(ix+2,self.Nx)):
                for j in range(max(iy-1,0),min(iy+2,self.Ny)):
                    for k in range(max(iz-1,0),min(iz+2,self.Nz)):
                        for a in self.grid[i][j][k]:
                            result.add(a)
            return result

        for i,a in enumerate(listOfAtoms):
            ix, iy, iz = self._find_grid_(a.x, a.y, a.z)
            self.neighborList[i] = allAtomsInNeighborsOfAGridBlock(ix,iy,iz)
            self.neighborList[i].remove(a)  # remove itself

        #self._checking_()

    def GetNeighborList(self,index):
        # returns a copy of the neighborlist. Just a copy, don't modify it
        return self.neighborList[index]

    def _checking_(self):
        from math import sqrt
        for i,a in enumerate(self.atomList):
            neighlist = self.GetNeighborList(i)
            output("Neighbor List of atom {} contains {} atoms.".format(i,len(neighlist)))
            for b in neighlist:
                dist = Distance(a,b)
                if dist > sqrt(3)*2*self.gridSize:
                    error("Something is wrong in _check_, distance = {}".format(dist))



    def _find_grid_(self,x,y,z):
        # return a list [ix,iy,iz] indicating the grid index of x, y, z
        # ix, iy, iz must be within [0, Nx], [0, Ny], [0, Nz]. Otherwise something is wrong.
        from math import floor
        ix = int(floor((x - self.minx) / self.gridSize))
        iy = int(floor((y - self.miny) / self.gridSize))
        iz = int(floor((z - self.minz) / self.gridSize))

        # Debugging, double-checking:
        # if ix < 0 or ix >= self.Nx or iy < 0 or iy >= self.Ny or iz < 0 or iz >= self.Nz:
        #     error("Something is wrong within _find_grid_()")

        return [ix,iy,iz]

    def DetectClashing(self,guestAtomList,minDist):
        # This function tests whether atoms in the guestAtomList clashes with any of the atoms that are used to
        # construct this neighbor list. minDist is the minimal distance that is considered to be a clash.
        # This function, if successful, returns a list of Bool with the same length as the guestAtomList, indicating
        # whether or not each atom in guestAtomList clashes with the host system. (True for clash, False for non-clash)
        # If a guestAtom is outside the box of the host system, ie, guestAtom.x not in [self.minx, self.maxx], etc, it's
        # considered to be a non-clash.
        result = [None for i in range(len(guestAtomList))]
        neighbors = None
        for index, guestAtom in enumerate(guestAtomList):
            result[index] = False
            ix,iy,iz = self._find_grid_(guestAtom.x,guestAtom.y,guestAtom.z)
            if ix<0 or ix>=self.Nx or iy<0 or iy>= self.Ny or iz<=0 or iz>=self.Nz:
                # guest atom outside the box, it's a non-clash
                continue
            neighbors = self.grid[ix][iy][iz]
            for hostAtom in neighbors:
                dist = Distance(hostAtom,guestAtom)
                if dist <= minDist:
                    result[index] = True
                    break

        return result






class Rule:
    ele1: str
    ele2: str
    low: float
    high: float
    type: str
    def __init__(self):
        pass

class BondRules(BondDetector):
    def __init__(self,globalCutoff = 5.0):
        self.rules = []
        self.globalCutoff = globalCutoff
        # Optional parameter global cutoff is used to build neighbor list as the grid size. Set this parameter to
        # be slightly larger than the maximum bond length in your system. If this value is too large, bond detection
        # may be very slow for large systems. If this value is too small, you may miss some bonds.



    def ParseFile(self,file_name):
        lines = None
        try:
            with open (file_name,"r") as file:
                lines = file.readlines()
        except:
            error("Can't open bond rule file <{}>".format(file_name),False)
            return False

        for i,line in enumerate(lines, start=1):
            l = line.strip(',')
            if (len(l)==0):
                continue
            rule = Bond()
            try:
                parts = l.split(',')
                rule.ele1 = parts[0]
                rule.ele2 = parts[1]
                rule.low = float(parts[2])
                rule.high = float(parts[3])
                rule.type = parts[4]
            except:
                error("Unexpected format @ line {} of bond rule file <{}>:\n{}".format(i,file_name,l),False)
                return False
            self.rules.append(rule)

        return True

    def Detect(self,molecularSystem,flushCurrentBonds):
        # It is required that the systemwideSerial of each atom must be unique.
        systemwideSerialToAtomMap = {}
        systemwideSerialToMoleculeMap = {}
        atomList = []
        tempBondList = []


        if flushCurrentBonds:
            molecularSystem.interMolecularBonds = []

        for m in molecularSystem.molecules:
            if flushCurrentBonds:
                m.bonds = []
            for a in m.atoms:
                if a.systemwideSerial == None or a.systemwideSerial in systemwideSerialToAtomMap:
                    error("In BondRule:Detect(), each atom's systemwideSerial must be unique", True)
                    return False
                systemwideSerialToAtomMap[a.systemwideSerial] = a
                systemwideSerialToMoleculeMap[a.systemwideSerial] = m
                atomList.append(a)

        # Performs a 2-body scan.
        # For system contains more than 10^4 atoms, O(N^2) is not acceptable. A neighbor list is needed here.
        N = len(atomList)
        if N > 10000:
            sys.stdout.write("Building NeighborList for system {} with grid size {} Å...".format(molecularSystem.name,self.globalCutoff))
        neighList = NeighborList(atomList,self.globalCutoff)
        if N > 10000:
            sys.stdout.write("Done.\n")
            sys.stdout.write("Scanning for Bonds:\n")

        for i in range(N):
            if N > 10000 and i%1000 == 0:
                ProgressBar(float(i)/N)

            a1 = atomList[i]
            neighbors = neighList.GetNeighborList(i)
            for a2 in neighbors:
                # To prevent the same bond from appearing twice, such as A-B and B-A, we require that the
                # systemwideSerial of a2 being greater than a1.
                # What are now comparing are two strings rather than integers. But it serves our purpose.
                if a1.systemwideSerial >= a2.systemwideSerial:
                    continue

                distance = Distance(a1,a2)
                newBond = None
                for r in self.rules:
                    if (a1.element != r.ele1 or a2.element != r.ele2) and (a1.element != r.ele2 or a2.element != r.ele1):
                        continue
                    if distance < r.low or distance > r.high:
                        continue
                    # Now we have a confirmed bond
                    newBond = Bond()
                    newBond.atom1 = a1.systemwideSerial
                    newBond.atom2 = a2.systemwideSerial
                    newBond.type = r.type
                    newBond.length = distance
                # Note that through the search of bond rules, multiple rules may be satisfied, but only the last satisfied
                # rule will be recorded. In other words, bond rules that appear later will override previous ones.
                if newBond != None:
                    tempBondList.append(newBond)
        if N > 10000:
            ProgressBar(1.0)
            output('')



        # Now we add the bonds to the system.
        # If old bonds are not removed, extra work is need to make sure that this are no duplicate bonds.
        existingBonds = set()
        if not flushCurrentBonds:
            # Going to construct a set to record which bonds are already present.
            # The 'existingBonds' is a set that records pairs of atoms that are bonded previously.
            # If A and B are bonded, then both 'A.#-B.#' and 'B.#-A.#' are stored in it, where # means atom's systemwideSerial
            for m in molecularSystem.molecules:
                for b in m.bonds:
                    a1 = None
                    a2 = None
                    for a in m.atoms:
                        if a.serial == b.atom1:
                            a1 = a
                        elif a.serial == b.atom2:
                            a2 = a
                        else:
                            pass
                    key1 = '{}-{}'.format(a1.systemwideSerial,a2.systemwideSerial)
                    key2 = '{}-{}'.format(a2.systemwideSerial,a1.systemwideSerial)
                    existingBonds.add(key1)
                    existingBonds.add(key2)
            for b in molecularSystem.interMolecularBonds:
                key1 = '{}-{}'.format(b.atom1,b.atom2)
                key2 = '{}-{}'.format(b.atom2,b.atom1)
                existingBonds.add(key1)
                existingBonds.add(key2)



        totalBonds = len(tempBondList)
        actuallyAddedBonds = 0

        for bond in tempBondList:
            if not flushCurrentBonds:
                key = '{}-{}'.format(bond.atom1, bond.atom2)
                if key in existingBonds:
                    continue

            a1 = systemwideSerialToAtomMap[bond.atom1]
            a2 = systemwideSerialToAtomMap[bond.atom2]
            m1 = systemwideSerialToMoleculeMap[bond.atom1]
            m2 = systemwideSerialToMoleculeMap[bond.atom2]
            if m1 == m2: # Intermolecular bond
                bond.atom1 = a1.serial
                bond.atom2 = a2.serial
                m1.bonds.append(bond)
            else: # Intramolecular bonds
                molecularSystem.interMolecularBonds.append(bond)

            actuallyAddedBonds += 1

        # For debugging:
        output('Totally {} bonds found, {} bonds are added.'.format(totalBonds,actuallyAddedBonds))

        return True

def DefaultBondRules(globalCutoff = 5.0):
    # This is a factory method that returns a BondRules object which reads the default bond rule file 'bondrules.csv'
    rules = BondRules(globalCutoff)
    import os
    default_route = '{}/{}'.format(os.path.dirname(__file__),"bondrules.csv")
    # the first part of the path is needed if case the script is called from other directories
    result = rules.ParseFile(default_route)
    return rules

def TestBondDetection():
    r = DefaultBondRules()
    ms = MolecularSystem()

    ms.Read(PDBFile(),'GenerateMoltemplateLTFileFromMOL2/000.Common.Files.Read.Only/water.and.ch4.cells/4158methane.pdb')


    ms2 = MolecularSystem()
    # ms2 = BreakupMoleculeRandomly(ms.molecules[0],10)
    # ms2.interMolecularBonds = []
    ms2 = ms

    # randomly delete some bonds:
    def RandDelBonds():
        for m in ms2.molecules:
            import random
            random.shuffle(m.bonds)
            for i in range(random.randint(0,len(m.bonds))):
                del(m.bonds[-1])
        random.shuffle(ms2.interMolecularBonds)
        for i in range(random.randint(0,len(ms2.interMolecularBonds))):
            del(ms2.interMolecularBonds[-1])

    for i in range(0):
        RandDelBonds()
        ms2.AutoDetectBonds(r)

    ms2.Summary()
    ms3 = MolecularSystem()
    ms3.molecules.append( ReduceSystemToOneMolecule(ms2))
    ms3.Summary()
    ms3.RenumberAtomSerials()
    output.setoutput(open('testcase/dump.mol2','w'))
    ms3.Write(MOL2File())

if __name__ == 'main':
    TestBondDetection()