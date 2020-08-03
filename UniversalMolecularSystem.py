#!/usr/bin/python
import sys
import math
class Bond:
    def __init__(self):
        self.atom1 = None
        self.atom2 = None
        self.low = 0.0
        self.high = 0.0
        self.order = None
        self.length = 0.0
    def Read(self,line):
        parts = line.split()
        if (len(parts) != 5):
            return False
        self.atom1 = parts[0]
        self.atom2 = parts[1]
        self.low = float(parts[2])
        self.high = float(parts[3])
        self.order = parts[4]
        return True
    def Satisfy(self,atom1,atom2,dist):
        if ( (atom1 == self.atom1 and atom2 == self.atom2) or 
             (atom1 == self.atom2 and atom2 == self.atom1) ):
            if float(dist) > self.low and float(dist) < self.high:
                return True
        return False

    def Show(self):
        print(str(self.atom1) + " " + str(self.atom2) + " " + str(self.low)+ " " +
              str(self.high) + " " + str(self.order)+" actual length = "+str(self.length))

class Atom:
    def __init__(self):
        self.name = None   # Its element
        self.x = self.y = self.z = 0.0
        self.internalName = None    # Its internal name as appears in mol2 files, like "H13" or "C2"
        self.type = None   # for example its hybridization type, like "C.3" or "P.2"
        self.charge = 0
        self.flexible = True
        self.layerInfo = ""
    def ShowAsXYZ(self):
        print(str(self.name) + " " +
              str(self.x) + " " +
              str(self.y) + " " +
              str(self.z) )
    def ShowAllFields(self):
        print(str(self.name) + " " +
              str(self.x) + " " +
              str(self.y) + " " +
              str(self.z) + " " +
              str(self.internalName) + " " +
              str(self.type) + " " +
              str(self.flexible)
              )

class Molecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
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
                    bond.order = order
                    self.bonds.append(bond)

    def ParseXYZFileLines(self,lines):
        for i in range(2,len(lines)):
            pars = lines[i].split()
            if (len(pars) < 4):
                return False
            a = Atom()
            a.name = pars[0]
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

    def ParseMol2FileLine(self,line):
        if len(line.strip()) == 0:
            return True
        pars = line.split()
        if len(pars) < 6:
            return False
        a = Atom()
        a.name = pars[5].split(".")[0]
        try:
            a.x = float(pars[2])
            a.y = float(pars[3])
            a.z = float(pars[4])
        except ValueError:
            return False
        a.internalName = pars[1]
        a.type = pars[5]
        self.atoms.append(a)
        return True


    def ShowAsXYZ(self):
        print("  "+str(len(self.atoms)))
        print("xyz2mo2.py")
        for i in self.atoms:
            i.ShowAsXYZ()

    def ShowAsMol2(self):
        print("@<TRIPOS>MOLECULE")
        print("mol")
        print(str(len(self.atoms))+" "+str(len(self.bonds)))
        print("SMALL")
        print("NO_CHARGES")
        print("")
        print("")
        print("@<TRIPOS>ATOM")
        for i in range(0,len(self.atoms)):
            print(str(i+1) + "  " + str(self.atoms[i].name) + " " +
                  str(self.atoms[i].x) + " " + 
                  str(self.atoms[i].y) + " " +
                  str(self.atoms[i].z) + " " +
                  str(self.atoms[i].name) )
        print("@<TRIPOS>BOND")
        for i in range(0,len(self.bonds)):
            b = self.bonds[i]
            print(str(i+1)+" "+str(b.atom1+1)+" "+str(b.atom2+1)+" "+str(b.order))

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

class Mol2File:  # Its main function "ParseFile()" acts like a factory for "Molecule"s
    def __init__(self):
        self.name = None
        self.molecules = []

    def ParseFile(self, file_name):

        from enum import Enum
        class ReadingStatus(Enum):
            Atom = 1
            Bond = 2
            Other = 3

        file = None
        try:
            file = open(file_name, "r")
        except IOError:
            sys.stderr.write("Can't open mol2 file <" + str(file_name) + "> \n")
            return False
        content = file.readlines()
        i = 0

        mol = None
        fileStatus = ReadingStatus.Other

        while i < len(content):
            line = content[i]
            if line.startswith("@<TRIPOS>MOLECULE"):
                if mol!= None:  # Found the 1st molecule ( == None) or found a new molecule ( != None)
                    self.molecules.append(mol)
                mol = Molecule()
                fileStatus = ReadingStatus.Other

            elif line.startswith("@<TRIPOS>ATOM"):
                fileStatus = ReadingStatus.Atom

            elif line.startswith("@<TRIPOS>BOND"):
                fileStatus = ReadingStatus.Bond

            elif line.startswith("@"):
                fileStatus = ReadingStatus.Other
            else:
                pass

            if fileStatus == ReadingStatus.Atom:
                mol.ParseMol2FileLine(line)
            elif fileStatus == ReadingStatus.Bond:
                b = Bond()
                pars = line.split()
                if(len(pars) == 4):
                    b.atom1 = int(pars[1])
                    b.atom2 = int(pars[2])
                    b.order = pars[3]
                    mol.bonds.append(b)
            else:
                pass

            i = i+1

        self.molecules.append(mol) # Append the last molecule

        if (len(self.molecules) == 0):
            sys.stderr.write("Error occurred while reading the mol2 file <" + str(file_name) + ">, no"
                             "molecular structure was read in.")
            return False
        file.close()
        return True

class GaussianLogFile:  # Its main function "ParseFile()" acts like a factory for "Molecule"s
    def __init__(self):
        self.name = None
        self.molecules = []

    def ParseFile(self,file_name):

        PeriodicTable="H He "\
        "Li Be B C N O F Ne "\
        "Na Mg Al Si P S Cl Ar "\
        "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "\
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "\
        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "\
        "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo"
        elements = PeriodicTable.split()

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
                a.name = elements[int(parts[1])-1]
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

def MainAsMol22XYZ(argc,argv):
    if argc != 2:
        sys.stderr.write("Usage : mol22xyz.py mol2file")
        exit()

    file = Mol2File()
    file.ParseFile(argv[1])
    for m in file.molecules:
        m.ShowAsXYZ()

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
                        that appeared eariler, if it controdicts with the earlier rules.
                        Other specifications include:
                        1. In each rule, the order of names of two atoms doesn't matter, meaning "C O 1.2 1.6 1" and "O C 1.2 1.6 1" are
                           exactly the same.
                        2. The 3rd parameter of this program, [bond_rule_file], is optional. There is a default bond rule file, named
                           "xyz2mol2.default.bond.rules", comes together with the xyz2mol2.py program and read by the program each time
                           it runs, to recognize some commonly encountered bonds in organic molecules. However a user defined [bond_rule_file]
                           can be supplied if some bonds are not recognized or special rules apply. Remember that rules in [bond_rule_file]
                           will override the default rules, if there are controdictions.
                        3. The 5th segment of each line in bond_rule_file defines the bond type, it can be any string that conforms with the
                           TRIPOS mol2 file format specification, such as 1, 2, 3, ar, am, etc. More information can be obtained from
                           www.tripos.com.
            """)
        exit()

    rules = BondRules()
    default_route = "/home/nfs/luowj/bin/xyz2mol2.default.bond.rules"
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

if __name__ == '__main__':
    sys.argv = ["","clean.xyz","INCAR","0 1 0 1","0 1","TFInfo","LayerInfo"]
    MainAsXYZ2Gaussian(len(sys.argv),sys.argv)
