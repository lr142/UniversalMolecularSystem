# Automatically detect chemical bonds based on element types and distances.
# It should be implemented in a more advanced way, for example through machine learning...
# I'll work on this later...

# The class BondRules and the function MainAsXYZ2Mol2 are from legacy code, just for a reference.
# The class Bond is an older version of the Bond Class, which should also be rewritten.


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
        i.ShowAsMol2()1000

class Bond:
    def __init__(self):
        self.atom1 = None
        self.atom2 = None
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