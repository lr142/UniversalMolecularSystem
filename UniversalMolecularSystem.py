#!/usr/bin/python
import sys
import math
from Utility import *

class Bond:
    def __init__(self):
        self.atom1 = None
        self.atom2 = None
        self.low = 0.0
        self.high = 0.0
        self.type = None
        self.length = 0.0
    def Read(self,line):
        parts = line.split()
        if (len(parts) != 5):
            return False
        self.atom1 = parts[0]
        self.atom2 = parts[1]
        self.low = float(parts[2])
        self.high = float(parts[3])
        self.type = parts[4]
        return True
    def Satisfy(self,atom1,atom2,dist):
        if ( (atom1 == self.atom1 and atom2 == self.atom2) or 
             (atom1 == self.atom2 and atom2 == self.atom1) ):
            if float(dist) > self.low and float(dist) < self.high:
                return True
        return False

    def Show(self):
        print(str(self.atom1) + " " + str(self.atom2) + " " + str(self.low) + " " +
              str(self.high) + " " + str(self.type) + " actual length = " + str(self.length))

class Atom:
    element: str
    x:float
    y:float
    z:float
    type:str
    charge:float
    flexible:bool
    serial:str
    layerInfo: str

    def __init__(self):
        self.element = None   # Its element, must be correctly capitalized like "C" or "Rh"; "C1", "c", "rh" are not acceptable.
        self.x = self.y = self.z = 0.0  # Cartesian or fractional coordinates
        self.name = None   # Its internal name as appears in mol2 files, like "H13" or "C2"
        self.type = None   # Its force field type or its hybridization type, like "C.3" or "P.2"
        self.charge = 0
        self.flexible = True
        self.serial = None   # serial number in some cases, like in PDB.
        self.layerInfo = None  # Used in QM/MM models, ie. high/mid/low layer in ONIOM methods
        # Not all properties are available or relevant in all cases. File readers/writers and users can add additional
        # properties to an Atom if necessary.
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
    name: str
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.name = None
        self.bondedTo = None
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

    def ParseXYZFileLines(self,lines):
        for i in range(2,len(lines)):
            pars = lines[i].split()
            if (len(pars) < 4):
                return False
            a = Atom()
            a.element = pars[0]
            try:
                a.x = float(pars[1])
                a.y = float(pars[2])
                a.z = float(pars[3])

                if len(pars) > 4:   # Additional keywords in .xyz file to specify selective dynamics
                    if pars[4].upper()[0] == 'F':
                        a.flexible = False

            except ValueError:
                return False

            self.atoms.append(a)
        return True

    def ShowAsXYZ(self):
        print("  "+str(len(self.atoms)))
        print("xyz2mo2.py")
        for i in self.atoms:
            i.ShowAsXYZ()

    def Check(self, renameAtoms = True):
        # If all checks are passed, return True, otherwise False
        # This is a consistency check to make sure that:
        # 1. Atom serial numbers start from 1 and are consecutive (Not mandatory), or at least unique (mandatory)
        serialToIndexMap = {}
        warnedAboutConsecutiveIssues = False
        for i,atom in enumerate(self.atoms):
            serial = atom.serial
            if i == 0 and atom.serial != '1':   # Not start from 1
                error("Atom serial in the molecule doesn't start from 1, but from {}".format(serial),False)
            if not warnedAboutConsecutiveIssues:
                try:    # Not consecutive or the serial are not even numbers. This warning shall be given only once
                    if i > 0 and int(serial) != int(self.atoms[i-1].serial) + 1:
                        error("Atom serials in the molecule are not consecutive from {} to {}".format(
                            self.atoms[i-1].serial,serial),False)
                        warnedAboutConsecutiveIssues = True
                except ValueError:
                    error("Atom serial is not a number: {} -> {}".format(self.atoms[i-1].serial,serial),False)
                    warnedAboutConsecutiveIssues = True

            if serial in serialToIndexMap:  # serial are not unique
                error("Atom serials in the molecule are not unique, found that serial {} appears again".format(serial),False)
                return False

            serialToIndexMap[serial] = i

        # 2. Atom names are unique (This may not be true in many cases), if renameAtoms is set to True, atoms will
        #       be renamed to Element+Serial
        alreadyUsedAtomNames = set([])
        warnedAboutAtomNamingUniquenessIssue = False
        for atom in self.atoms:
            if (not warnedAboutAtomNamingUniquenessIssue) and (atom.name in alreadyUsedAtomNames):
                error("Atom names in the molecule are not unique, found"
                "that name {} appears again".format(atom.name),False)
                error("Atoms {} be renamed. This message will not appear again.".format(
                    "WILL" if renameAtoms else "WILL NOT"),False)
                warnedAboutAtomNamingUniquenessIssue = True

            alreadyUsedAtomNames.add(atom.name)

        if warnedAboutAtomNamingUniquenessIssue and renameAtoms:
            for atom in self.atoms:
                atom.name = "{}{}".format(atom.name,atom.serial)

        # 3. All bonds are correctly specified, i.e no bonds pointing to atoms that are not exist.
        #       In this step, a auxiliary data structure 'bondedTo' is created, which is an array of maps
        #       recording all other atoms (in their indexes) bonded to the current one
        self.bondedTo = [set() for i in range(len(self.atoms))]
        for b in self.bonds:
            fromSerial = b.atom1
            toSerial = b.atom2
            if fromSerial in serialToIndexMap and toSerial in serialToIndexMap:
                fromIndex = serialToIndexMap[fromSerial]
                toIndex = serialToIndexMap[toSerial]
                self.bondedTo[fromIndex].add(toIndex)
                self.bondedTo[toIndex].add(fromIndex)
            else:
                error("Dangling bond found in molecule: {} - {} with type {}".format(b.atom1,b.atom2,b.type),False)
                self.bondedTo = None  # Destroys the bondedMap and quit
                return False

        return True





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

class MolecularSystem:
    # A molecular system is a collection of molecules and, in some cases, associated information including boundary
    # conditions, force field descriptions, overall thermodynamic properties, etc.
    # It can also represent the trajectory of a MD simulation if the system contains only 1 molecule, for example in
    # analyzing the results of a QM run. If the system contains multiple molecules, a collection of MolecularSystems
    # shall be needed to represent a trajectory.
    def __init__(self):
        self.molecules = []
        self.boundary = None;    # boundary should be a 3x3 matrix for an orthogonal system
    def Read(self,molecularFile,filename):
        molecularFile.Read(self,filename)
    def Write(self,molecularFile):
        molecularFile.Write(self)
    # Read and Write operations are delegated to concrete MolecularFile classes.
    def Summary(self):
        molCount = len(self.molecules)
        atomCount = 0
        bondCount = 0
        for m in self.molecules:
            atomCount += len(m.atoms)
            bondCount += len(m.bonds)
        output("Molecular System has {} molecules, {} atoms, and {} bonds".format(molCount,atomCount,bondCount))



class XYZFile:  # Its main function "ParseFile()" acts like a factory for "Molecule"s
    def __init__(self):
        self.name = None
        self.molecules = []

    def ParseFile(self,file_name):

        file = None
        try:
            file = open(file_name,"r")
        except IOError:
            sys.stderr.write("Can't open xyz file <"+str(file_name)+ "> \n")
            return False

        content = file.readlines()
        
        i = 0
        while i < len(content):
            pars = content[i].strip().split()
            if(len(pars)==0):
                i = i+1
                continue
            
            atom_count = None
            try:
                atom_count = int(pars[0])
            except ValueError:
                sys.stderr.write("Error occurred while reading line " +
                                 str(i+1) + " from file <" + str(file_name) + ">, an atom count is expected.\n")
                return False
                
            end_line_no = i+2+atom_count

            if(atom_count == 0 or end_line_no > len(content)):
                sys.stderr.write("Error occurred while reading the molecule starting from line " +
                                 str(i+1) + " from file <" + str(file_name) + ">, check the atom count.\n")
                return False

            
            m = Molecule()
            result = m.ParseXYZFileLines(content[i:end_line_no])

            if (result):
                self.molecules.append(m)
            else:
                sys.stderr.write("Error occurred while reading the molecule starting from line " +
                                 str(i+1) + " from file <" + str(file_name) + ">\n")
            i = end_line_no

        if(len(self.molecules) == 0):
            sys.stderr.write("Error occurred while reading the xyz file <" + str(file_name) +">, no "
                             "molecular structure was read in.")
            return False


        file.close()
        return True

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

class BondRules:
    def __init__(self):
        self.rules = []

    def ParseFile(self,file_name):
        file = None
        try:
            file = open(file_name,"r")
        except IOError:
            sys.stderr.write("Can't open bond rule file <"+ str(file_name)+">/n")
            return False

        lines = file.readlines()

        for i in range(0,len(lines)):
            l = lines[i].strip()
            if (len(l)==0):
                continue
            rule = Bond()
            result = rule.Read(l)
            if (not result):
                sys.stderr.write("Error in parsing <"+str(file_name+"> at line "+str(i+1)+" :\n"))
                sys.stderr.write(lines[i]+"\n")
                sys.stderr.write("Format error.\n")
                return False
            self.rules.append(rule)
        file.close()
        return True

    def Show(self):
        for i in self.rules:
            i.Show()
    def CheckRules(self,a1,a2,dist):
        for i in reversed(self.rules):
            if ( i.Satisfy(a1, a2, dist)):
                return i.order
        return None

def MainAsXYZ2Mol2(argc,argv):
    if argc < 2 or argc > 3:
        sys.stderr.write(
            """
                        Usage : xyz2mol2.py xyzfile [bond_rule_file]
                        The [bond_rule_file] is a txt files which defines bond rules, they can include multiple lines, each line is like below:
                        AtomType1 AtomType2 MinBondLength MaxBondLength BondType
                        For example in the following lines
                        C O 1.2 1.6 1
                        C O 1.1 1.2 ar
                        The above two lines ask the program to judge any two C, O atoms that are within 1.2 Angstrom(A) but at least 1.1 A
                        apart form a aromatic bond, any two C, O atoms that are at least 1.2 A apart but within 1.6 A from each other form
                        a single bond. Rules are checked in reverse order, which means that a rule appears later in the file overrides rules
                        that appeared earlier, if it contradicts with the earlier rules.
                        Other specifications include:
                        1. In each rule, the order of names of two atoms doesn't matter, meaning "C O 1.2 1.6 1" and "O C 1.2 1.6 1" are
                           exactly the same.
                        2. The 3rd parameter of this program, [bond_rule_file], is optional. There is a default bond rule file, named
                           "xyz2mol2.default.bond.rules", comes together with the xyz2mol2.py program and read by the program each time
                           it runs, to recognize some commonly encountered bonds in organic molecules. However a user defined [bond_rule_file]
                           can be supplied if some bonds are not recognized or special rules apply. Remember that rules in [bond_rule_file]
                           will override the default rules, if there are contradictions.
                        3. The 5th segment of each line in bond_rule_file defines the bond type, it can be any string that conforms with the
                           TRIPOS mol2 file format specification, such as 1, 2, 3, ar, am, etc. More information can be obtained from
                           www.tripos.com.
            """)
        exit()

    rules = BondRules()
    default_route = "bond.rules"
    result = rules.ParseFile(default_route)
    if (not result):
        print("Cannot find the default bond rules file, hope you know what you are doing...")
        print("The default bond rules file should be located at:")
        print(default_route)
        pass


    if (argc == 3):
        result = rules.ParseFile(argv[2])
        if (not result):
            exit()

    #    rules.Show()

    xyzfile = XYZFile()
    xyzfile.ParseFile(argv[1])

    for i in xyzfile.molecules:
        i.FindBonds(rules)

    for i in xyzfile.molecules:
        i.ShowAsMol2()

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


