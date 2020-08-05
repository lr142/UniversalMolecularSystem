from UniversalMolecularSystem import *
from Utility import *

class MOL2File(MolecularFile):
    def __init__(self):
        pass

    def Read(self,molecularSystem,filename):

        from enum import Enum
        class ReadingStatus(Enum):
            Atom = 1
            Bond = 2
            Other = 3
            CRYSIN = 4
            Molecule = 5

        file = None
        try:
            file = open(filename, "r")
        except IOError:
            error("Can't open mol2 file <" + str(self.filename) + "> \n")
            return False
        content = file.readlines()
        file.close()
        i = 0

        molecularSystem.molecules = []
        mol = None
        fileStatus = ReadingStatus.Other

        while i < len(content):
            line = content[i]

            if line.startswith("@"):
                # Special Character that means the start of a new section
                if line.startswith("@<TRIPOS>MOLECULE"):
                    if mol != None:  # Found the 1st molecule ( == None) or found a new molecule ( != None)
                        molecularSystem.molecules.append(mol)
                    mol = Molecule()
                    fileStatus = ReadingStatus.Molecule

                elif line.startswith("@<TRIPOS>ATOM"):
                    fileStatus = ReadingStatus.Atom

                elif line.startswith("@<TRIPOS>BOND"):
                    fileStatus = ReadingStatus.Bond
                elif line.startswith("@<TRIPOS>CRYSIN"):
                    # This is the default way Materials Studio is used to write the periodic boundary condition.
                    fileStatus = ReadingStatus.CRYSIN
                elif line.startswith("@"):
                    fileStatus = ReadingStatus.Other
                else:
                    pass
                    # more keywords and more features can be read from the MOL2 file. If necessary, extend the program
                    # from the above enumeration statements and add respective actions below.

                i = i+1
                continue   # jump to the next line

            # if reaches here, the this line doesn't start with a "@"
            if fileStatus == ReadingStatus.Molecule:
                molName = line.strip().strip("\"")
                mol.name = molName
                fileStatus = ReadingStatus.Other

            elif fileStatus == ReadingStatus.Atom:
                a = self.ParseAtomLine(line)
                if a != None:
                    mol.atoms.append(a)
            elif fileStatus == ReadingStatus.Bond:
                b = Bond()
                pars = line.split()
                if (len(pars) == 4):
                    b.atom1 = pars[1]
                    b.atom2 = pars[2]
                    b.type = pars[3]
                    mol.bonds.append(b)
            elif fileStatus == ReadingStatus.CRYSIN:
                #Format: lx ly lz A B C space_group setting
                # Here I'm assuming it's an orthogonal cell
                pars = line.split()
                try:
                    lx = float(pars[0])
                    ly = float(pars[1])
                    lz = float(pars[2])
                except:
                    error("Unexpect format in <" + str(self.filename) + "> \n{}".format(line))
                    return False
                molecularSystem.boundary = [ [lx,0,0], [0,ly,0], [0,0,lz]]
            else:
                pass

            i = i + 1

        molecularSystem.molecules.append(mol)  # Append the last molecule

        if (len(molecularSystem.molecules) == 0):
            error("Error occurred while reading the mol2 file <" + str(self.filename) + ">, no"
                  "molecular structure was read in.",False)
            return False
        file.close()
        return True

    def ParseAtomLine(self,line):
        if len(line.strip()) == 0:
            return None
        pars = line.split()

        a = Atom()
        try:
            # Format:
            # atom_id atom_name x y z atom_type [subst_id
            # [subst_name [charge [status_bit]]
            # where "subst**" are information we don't care.
            # Get the charge from pars[8] if it has such information
            a.serial = pars[0]  # serial is read as a str
            a.name = pars[1]
            a.x = float(pars[2])
            a.y = float(pars[3])
            a.z = float(pars[4])
            a.name = pars[1]
            a.type = pars[5]
            a.element = a.type.split(".")[0]
            a.charge = float(pars[8]) if len(pars) > 8 else 0.0

        except:
            error("Unexpected format in mol2 file <" + str(self.filename) + ">, "
            "while reading the line:\n{}".format(line))

        return a

    def Write(self,molecularSystem):
        def WriteAMolecule(mol):
            output("@<TRIPOS>MOLECULE")
            output(mol.name if mol.name != None else "mol")
            output(str(len(mol.atoms)) + " " + str(len(mol.bonds)))
            output("SMALL")
            output("USER_CHARGES")
            output("")
            output("")
            output("@<TRIPOS>ATOM")
            for a in mol.atoms:
                output("{} {} {} {} {} {} {} {} {}".format(
                    a.serial,a.name,a.x,a.y,a.z,a.type,0.0,"****",a.charge)
                )
            output("@<TRIPOS>BOND")
            for i,b in enumerate(mol.bonds,start=1):
                print(str(i) + " " + str(b.atom1) + " " + str(b.atom2) + " " + str(b.type))

        for m in molecularSystem.molecules:
            WriteAMolecule(m)

        if molecularSystem.boundary != None:
            output("@<TRIPOS>CRYSIN")
            output("{} {} {} 90 90 90 1 1".format(
                molecularSystem.boundary[0][0],
                molecularSystem.boundary[1][1],
                molecularSystem.boundary[2][2]
            ))




def TestMOL2File(filename):

    s = MolecularSystem()

    s.Read(MOL2File(),filename)
    result = s.molecules[0].Check()
    output('Check {}'.format("passed" if result else "NOT PASSED!"))

    # file = open("dump.mol2",'w')
    # output.setoutput(file)
    # s.Write(MOL2File())

    s.Summary()
