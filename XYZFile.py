from UniversalMolecularSystem import *

from MOL2File import *

class XYZFile(MolecularFile):
    def __init__(self):
        pass

    def Read(self,molecularSystem,filename):
        # XYZ file usually contains only element and x,y,and z coordinates.
        # In our system we shall extend the XYZ format a little bit, to allow an trailing flag "T" or "F" to
        # denote whether an atom is flexible in the QM or MD simulation.
        # Note that an XYZ file can contain more than 1 molecule. It's also optional to supply a molecule name
        # in the comment line ( the 2nd line, or the line below each atom count if multiple molecules are present)

        file = None
        try:
            file = open(filename, "r")
        except IOError:
            error("Can't open xyz file <{}>".format(filename))
            return False
        content = file.readlines()
        file.close()

        from enum import Enum
        class ReadingStatus(Enum):
            Molecule = 1
            Atom = 2
            Other = 3

        status = ReadingStatus.Other

        i = 0
        atomCount = None
        atomIndex = None
        mol = None
        molecularSystem.molecules = []

        while i < len(content):
            line = content[i].strip()
            if status == ReadingStatus.Other:
                if len(line) == 0:   # Skip blank lines. This happens near the end of the file.
                    i += 1
                    continue
                # Assert that an integer (atom count) should be read
                status = ReadingStatus.Molecule
                try:
                    atomCount = int(line)
                    if atomCount <=0:
                        raise Exception()# Make sure that the atom count is not zero, otherwise it throws an exception
                except:
                    error("Unexpected file format in xyz file <{}> at the following line."
                          " An atom count (>0) is expected.\n{}".format(filename,line))
                    return False
                mol = Molecule()
                atomIndex = 0
                mol.atoms = [Atom() for i in range(atomCount)]
            elif status == ReadingStatus.Molecule:
                # Read the comment line. If it contains some information, regard it as the molecule name
                if len(line) > 0:
                    mol.name = line
                status = ReadingStatus.Atom
            elif status == ReadingStatus.Atom:
                # Read in the atom information
                pars = line.split()
                atom = mol.atoms[atomIndex]  # This creates a link to the atom
                try:
                    atom.name = atom.element = pars[0]
                    atom.x = float(pars[1])
                    atom.y = float(pars[2])
                    atom.z = float(pars[3])
                    if len(pars) > 4:
                        atom.flexible = True if pars[4] == 'T' else False
                    atom.serial = "{}".format(atomIndex+1)
                except:
                    error("Unexpected file format in xyz file <{}> at the following line.\n{}".format(filename,line))
                    return False
                atomIndex += 1
                if atomIndex == atomCount:   # Read all atoms of a molecule
                    molecularSystem.molecules.append(mol)
                    status = ReadingStatus.Other
            else:
                pass

            i = i+1
            continue

        result = None
        # After the loop, make sure that there are not missing atoms
        if status != ReadingStatus.Other:
            error("Unexpected file format in xyz file <{}> at the end."
            " {} atoms are claimed but only {} are provided.".format(filename,atomCount,atomIndex))
            result = False
        else:
            result = True

        molecularSystem.RenumberAtomSerials()
        return result

    def Write(self,molecularSystem):
        # Writing out an XYZ file is straightforward
        for mol in molecularSystem.molecules:
            output("{}".format(len(mol.atoms)))
            output("{}".format(mol.name if mol.name != None else ""))
            for atom in mol.atoms:
                if atom.flexible == None:
                    output("{} {} {} {}".format(atom.element,atom.x,atom.y,atom.z,""))
                else:
                    output("{} {} {} {} {}".format(
                        atom.element, atom.x, atom.y, atom.z, "T" if atom.flexible else "F"))


def TestXYZFile(filename):

    s = MolecularSystem()

    s.Read(XYZFile(),filename)
    result = s.molecules[0].CheckConsistency()
    output('Check {}'.format("passed" if result else "NOT PASSED!"))

    s.Write(XYZFile())

    # file = open("dump.mol2",'w')
    # output.setoutput(file)

    s.Write(MOL2File())

    s.Summary()