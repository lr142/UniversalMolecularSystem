#!/usr/bin/python
import sys
import math
from Utility import *

class Bond:
    atom1: str
    atom2: str
    type: str
    length: float
    def __init__(self):
        self.atom1 = None
        self.atom2 = None
        self.type = None
        self.length = 0.0

    def Copy(self):
        import copy
        return copy.copy(self)

class Atom:
    element: str
    x:float
    y:float
    z:float
    type:str
    charge:float
    flexible:bool
    serial:str
    systemwideSerial: str
    layerInfo: str
    parent: str

    # Note: If complex data structures just as lists or references to other objects are included in Atom in future
    # versions, remember that Atom.Copy() is used to duplicate an atom, which uses copy.copy() internally to
    # 'shallow' copy atoms. Make sure such operations are safe
    # and reasonable. For example, if an atom contains a link to a third object, keep in mind that the external object
    # will not be copied. This note applies also to 'Bond'.
    # Since both Atom and Bond are considered to be 'elementary' classes, try not to include references or complex stuffs in them.

    def __init__(self):
        self.element = None   # Its element, must be correctly capitalized like "C" or "Rh"; "C1", "c", "rh" are not acceptable.
        self.x = self.y = self.z = 0.0  # Cartesian or fractional coordinates
        self.name = None   # Its internal name as appears in mol2 files, like "H13" or "C2". No need to be unique.
        self.type = None   # Its force field type or its hybridization type, like "C.3" or "P.2"
        self.charge = 0.0
        self.flexible = None
        self.serial = None     # serial number in a molecular, usually a number starting from 1 but can be any string
                               # serial number must be unique within a molecule
        self.systemwideSerial = None  # unique serial number in the molecular system. Must be present if there are
                               # intermolecular bonds (including hydrogen bonds) within the system.
        self.layerInfo = None  # Used in QM/MM models, ie. high/mid/low layer in ONIOM methods
        self.parent = None     # Used in some cases to record its parent, for example the molecule or residue it belongs to

        # Not all properties are available or relevant in all cases. File readers/writers and users can add additional
        # properties to an Atom if necessary.

    def Copy(self):
        # Returns a copy of the current Atom. It uses 'shallow' copy. See the notes at the beginning of this class
        import copy
        return copy.copy(self)

    def Translate(self,dx,dy,dz):  # Self-evident
        self.x += dx
        self.y += dy
        self.z += dz

    def ShowAsXYZ(self):
        print(str(self.element) + " " +
              str(self.x) + " " +
              str(self.y) + " " +
              str(self.z))
    def ShowAllFields(self):
        print(str(self.serial) + " " +
              str(self.element) + " " +
              str(self.x) + " " +
              str(self.y) + " " +
              str(self.z) + " " +
              str(self.name) + " " +
              str(self.type) + " " +
              str(self.flexible) + " " +
              str(self.charge)
              )

class Molecule:
    atoms: [Atom]
    bonds: [Bond]
    name: str
    bondedTo: [map]
    serial: str
    type: str

    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.name = None  # Names are user-defined identifiers, e.g. it can be type+serial in a protein. No need to be unique
        self.serial = None # Such as its serial in a peptide chain, must be unique within a MolecularSystem
        self.type = None # Such as the 3-Letter code residue name GLY
    def FindBonds(self,rules):
        for i in range(0,len(self.atoms)):
            #sys.stderr.write("Atom No. {}\n".format(i))
            for j in range(i+1,len(self.atoms)):
                a1 = self.atoms[i]
                a2 = self.atoms[j]
                length = math.sqrt (  (a1.x-a2.x)**2 +
                                 (a1.y-a2.y)**2 +
                                 (a1.z-a2.z)**2)
                order = rules.CheckRules(a1.name,a2.name,length)
                if ( order != None ):
                    bond = Bond()
                    bond.atom1 = int(i)
                    bond.atom2 = int(j)
                    bond.length = length
                    bond.type = order
                    self.bonds.append(bond)

    def BondedMap(self):
        # An auxiliary function, returns a list of sets (called bondedTo), for example
        # bondedTo[0] records the indexes (not serials) of atoms that are connected to atom with index (not serial) 0
        # Returns None if there are dangling bonds.
        bondedMap = [set() for i in range(len(self.atoms))]

        serialToIndexMap = {}
        for i,atom in enumerate(self.atoms):
            serialToIndexMap[atom.serial] = i

        for b in self.bonds:
            fromSerial = b.atom1
            toSerial = b.atom2
            if fromSerial in serialToIndexMap and toSerial in serialToIndexMap:
                fromIndex = serialToIndexMap[fromSerial]
                toIndex = serialToIndexMap[toSerial]
                bondedMap[fromIndex].add(toIndex)
                bondedMap[toIndex].add(fromIndex)
            else:
                return None

        return bondedMap

    def CheckConsistency(self):
        # If all checks are passed, return True, otherwise False
        # This is a consistency check to make sure that:
        # 1. Atom serial numbers must be unique
        serialToIndexMap = {}
        for i,atom in enumerate(self.atoms):
            serial = atom.serial
            if serial in serialToIndexMap:  # serial are not unique
                error("Atom serials in the molecule are not unique, found that serial {} appears again".format(serial),False)
                return False

            serialToIndexMap[serial] = i

        # 2. No dangling bonds, i.e no bonds pointing to atoms that are not exist.
        #    Used the auxiliary function 'BondedMap'
        result = self.BondedMap()
        if result == None:
            error("Dangling bond found in molecule: {} - {} with type {}".format(b.atom1,b.atom2,b.type),False)
            self.bondedTo = None  # Destroys the bondedMap and quit
            return False

        return True

    def Copy(self):
        import copy
        newMol = copy.copy(self)   # A 'Shallow' copy that copies only references to its string and integer properties
        newMol.atoms = []
        newMol.bonds = []
        for a in self.atoms:
            newMol.atoms.append(a.Copy())
        for b in self.bonds:
            newMol.bonds.append(b.Copy())
        return newMol


    def Summary(self):
        output("Molecule Serial: {}, Name: {}, Type: {}, with {} atoms, and {} bonds".format(
            self.serial,self.name,self.type,len(self.atoms),len(self.bonds)
        ))


class MolecularFile:
    # This is an abstract class that represents a molecular system description file and the methods to Read/Write
    # the file. For each supported file type, such as .xyz, .mol2, .pdb, etc., a separate concrete class shall be
    # created which implements the Read() and Write() virtual functions.
    def __init__(self):  # Virtual function
        pass
    def Read(self,molecularSystem,filename):   # Virtual function
        # Remember that any "Read" operation (even if not successful) will reset the molecular system
        # In cases where there is a need to read a system from multiple files, read them into separate
        # MolecularSystems and merge them.
        pass
    def Write(self,molecularSystem):   # Virtual function
        pass
class BondDetector:
    # This is an abstract class that has a sole function called Detect, which recognize and generates bonds within
    # a molecular system. In 'BondDetection.py' we implemented a clumsy detector based solely on atom distances and element
    # types. In the future we can write fancier detectors like one based on machine learning.
    def Detect(self,molecularSystem,flushCurrentBonds):
        pass

class MolecularSystem:
    # A molecular system is a collection of molecules and, in some cases, associated information including boundary
    # conditions, force field descriptions, overall thermodynamic properties, etc.
    # It can also represent the trajectory of a MD simulation if the system contains only 1 molecule, for example in
    # analyzing the results of a QM run. If the system contains multiple molecules, a collection of MolecularSystems
    # shall be needed to represent a trajectory.
    molecules: [Molecule]
    boundary: []
    interMolecularBonds: [Bond]
    name: str
    def __init__(self,name=None):
        self.molecules = []
        self.boundary = None;    # boundary should be a 3x3 matrix for an orthogonal system
        self.interMolecularBonds = []
        self.name = name;
    def Read(self,molecularFile,filename):
        molecularFile.Read(self,filename)
        # perform a consistency check, and create necessary auxiliary data structures
        for m in self.molecules:
            if not m.CheckConsistency():
                return False

    def Write(self,molecularFile):
        molecularFile.Write(self)

    def Copy(self):
        import copy
        newMS = copy.copy(self)
        newMS.molecules = []
        newMS.interMolecularBonds = []
        for m in self.molecules:
            newMS.molecules.append(m.Copy())
        for b in self.interMolecularBonds:
            newMS.interMolecularBonds.append(b.Copy())
        newMS.boundary = copy.deepcopy(self.boundary)
        return newMS

    def AtomCount(self):
        atomCount = 0
        for m in self.molecules:
            atomCount += len(m.atoms)
        return atomCount

    def RenumberAtomSerials(self, startingSystemwideSerial = 1):
        # renumber atom serials and systemwideSerials. This is required in most cases where the atoms must have a
        # consecutive order so that other programs can read it.
        oldToNewSerialMap = [ {} for i in range(len(self.molecules)) ]
        oldToNewSystemwideSerialMap = {}

        counterInSystem = startingSystemwideSerial
        for i,m in enumerate(self.molecules):
            counterInMolecule = 1
            for a in m.atoms:
                newSerial = '{}'.format(counterInMolecule)
                newSystemwideSerial = '{}'.format(counterInSystem)
                oldToNewSerialMap[i][a.serial] = newSerial
                if a.systemwideSerial != None:
                    oldToNewSystemwideSerialMap[a.systemwideSerial] = newSystemwideSerial
                a.serial = newSerial
                a.systemwideSerial = newSystemwideSerial
                counterInMolecule += 1
                counterInSystem += 1

        for i,m in enumerate(self.molecules):
            for b in m.bonds:
                b.atom1 = oldToNewSerialMap[i][b.atom1]
                b.atom2 = oldToNewSerialMap[i][b.atom2]

        for b in self.interMolecularBonds:
            b.atom1 = oldToNewSystemwideSerialMap[b.atom1]
            b.atom2 = oldToNewSystemwideSerialMap[b.atom2]

    def AutoDetectBonds(self, autoBondDetector, flushCurrentBonds = False):
        autoBondDetector.Detect(self, flushCurrentBonds)

    def Translate(self,dx,dy,dz):
        # Molecule.Translate is not defined. Seems unnecessary.
        for m in self.molecules:
            for a in m.atoms:
                a.Translate(dx,dy,dz)

    def FractionalToCartesianCoordinates(self):
        # Now support only orthogonal systems.
        A = self.boundary[0][0]
        B = self.boundary[1][1]
        C = self.boundary[2][2]
        for m in self.molecules:
            for a in m.atoms:
                a.x *= A
                a.y *= B
                a.z *= C

    # Read and Write operations are delegated to concrete MolecularFile classes.
    def Summary(self):
        molCount = len(self.molecules)
        bondCount = 0
        for m in self.molecules:
            bondCount += len(m.bonds)
            #m.Summary()
        output("MolecularSystem {} has {} molecules, {} atoms, {} intra-molecular bonds, and {} inter-molecular bonds".format(
            self.name,molCount,self.AtomCount(),bondCount,len(self.interMolecularBonds)))


# Legacy Codes to be rewritten
class GaussianLogFile:  # Its main function "ParseFile()" acts like a factory for "Molecule"s
    def __init__(self):
        self.name = None
        self.molecules = []

    def ParseFile(self,file_name):

        periodicTable = Utility.PeriodicTable()

        file = None
        try:
            file = open(file_name, "r")
        except IOError:
            sys.stderr.write("Can't open Gaussian Log file <" + str(file_name) + "> \n")
            return False
        content = file.readlines()

        from enum import Enum
        class ReadingStatus(Enum):
            AtomList = 0
            Atom = 1
            Other = 2

        mol = None
        fileStatus = ReadingStatus.Other
        mol = None
        i = 0
        while i < len(content):
            line = content[i]
            if line.strip() == "Input orientation:":
                fileStatus = ReadingStatus.Atom
                mol = Molecule()
                i = i + 5
                continue
            elif line.strip().startswith("-----------") and fileStatus == ReadingStatus.Atom:
                self.molecules.append(mol)
                fileStatus = ReadingStatus.Other
                i = i + 1
                continue
            else:
                i = i + 1

            if fileStatus == ReadingStatus.Atom:
                parts = line.split()
                x = float(parts[3])
                y = float(parts[4])
                z = float(parts[5])
                a = Atom()
                a.element = periodicTable.AtomicNumberToElement(int(parts[1]))
                a.x = x
                a.y = y
                a.z = z
                mol.atoms.append(a)

        if (len(self.molecules) == 0):
            sys.stderr.write("Error occurred while reading the Gaussian log file <" + str(file_name) + ">, no"
                             "molecular structure was read in.")
            return False
        file.close()
        return True

def MainAsXYZ2Gaussian(argc,argv):

    if (argc < 3 or argc > 7):
        sys.stderr.write("Usage : xyz2gaussianmol2.py xyzfile GaussianRouteFile [CHARGE] [SPIN_MULTIPLICITY]"
                         " [TFInfo] [LayerInfoForOniomModel]")
        exit()

    charge = 0
    spin_multiplicity = 1

    if argc > 3:
        charge = argv[3]
        spin_multiplicity = argv[4]

    xyzfile = XYZFile()
    xyzfile.ParseFile(argv[1])

    if argc >= 6:
        # Read the TFInfo file
        tfinfofile = open(sys.argv[5])
        for i in tfinfofile:
            l = i.strip()
            atomindex = int(l.split()[0]) - 1
            flexible = False if l.split()[1][0] == "F" else True
            # The TF flags should be applied to each molecule, for example IS, FS, and TS, equivalently.
            for m in xyzfile.molecules:
                m.atoms[atomindex].flexible = flexible

    if argc >= 7:
        # Read the Oniom Layer Info file
        layerinfofile = open(sys.argv[6])
        for i in layerinfofile:
            parts = i.strip().split()
            atomindex = int(parts[0])-1
            del(parts[0])
            content = " ".join(parts)
            # The Layer Info flags should be applied to each molecule, for example IS, FS, and TS, equivalently.
            for m in xyzfile.molecules:
                m.atoms[atomindex].layerInfo = content

    route = open(sys.argv[2], "r")
    for i in route:
        l = i.strip()
        if len(l) > 0:
            print(l)
    print("")

    index = 1
    for mol in xyzfile.molecules:
        print("Molecule {}".format(index))
        index += 1
        print("")
        print("{} {}".format(charge, spin_multiplicity))
        allAtomsFlexible = True
        for a in mol.atoms:
            if a.flexible == False:
                allAtomsFlexible = False
        for a in mol.atoms:
            if allAtomsFlexible:
                a.ShowAsXYZ()
            else:
                print("{}  {}   {}   {}   {}{}".format(
                    a.name,
                    0 if a.flexible else -1,
                    a.x,
                    a.y,
                    a.z,
                    "" if a.layerInfo == "" else " "+a.layerInfo
                ))
        print("")

def MainAsGaussianLog2XYZ(argc,argv):
    if argc != 2:
        sys.stderr.write("Usage : gaussianlog2xyz.py gaussian_log_file")
        exit()

    file = GaussianLogFile()
    file.ParseFile(argv[1])
    for m in file.molecules:
        m.ShowAsXYZ()


