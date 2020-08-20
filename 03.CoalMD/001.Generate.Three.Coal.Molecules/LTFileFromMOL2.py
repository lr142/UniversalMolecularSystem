import sys

sys.path.append("../..")

from Utility import *
from UniversalMolecularSystem import *
from MOL2File import *

class MoltemplateLTFile(MolecularFile):
    def __init__(self):
        self.writeBondInfo = True
        pass

    def AssignForcefieldType(self,molecule,force_field_type_file):
        with open(force_field_type_file,"r") as file:
            lines = file.readlines()

        serialToAtomMap = {}
        for atom in molecule.atoms:
            serialToAtomMap[atom.serial] = atom   # saves links to corresponding atoms
            atom.name = "{}{}".format(atom.element,atom.serial)  # need to rename the atoms into uniqueness in a force field file
        for l in lines:
            pars = l.strip().split(',')
            serial = pars[0]
            if serial == "":  # This csv file may contain irrelevant info. Just ignore them.
                return
            type = pars[1]
            serialToAtomMap[serial].type = type

    def Write(self,molSystem):
        mol = molSystem.molecules[0]
        # Assuming there is a single molecule in the system
        output('import "oplsaa.lt"')
        output('# Automatically generated from python script with user-designated force field type for each atom')
        output('{} inherits OPLSAA{{'.format(mol.name))
        output('  write(\'Data Atoms\'){')
        for a in mol.atoms:
            output('    $atom:{}   $mol:.   @atom:{}   0.0    {}   {}   {}'.format(
                a.name, a.type, a.x, a.y, a.z
            ))
        output('  }')
        if self.writeBondInfo:
            output('  write(\'Data Bond List\'){')

            serialToAtomMap = {}
            for atom in mol.atoms:
                serialToAtomMap[atom.serial] = atom   # saves links to corresponding atoms

            for i,b in enumerate(mol.bonds):
                a1 = serialToAtomMap[b.atom1]
                a2 = serialToAtomMap[b.atom2]
                output('    $bond:b{}   $atom:{}   $atom:{}'.format(i+1,a1.name,a2.name))
            output('  }')

        output('}')

file_names = [["Given1.mol2","Given1.csv","given1.lt","given1.nb.lt"],
            ["Given2.mol2","Given2.csv","given2.lt","given2.nb.lt"],
            ["FuchsSandoff.mol2","FuchsSandoff.csv","fuchssandoff.lt","fuchssandoff.nb.lt"],
              ["c4rings.mol2","c4rings.csv","c4rings.lt","c4rings.nb.lt"]

         ]

for i in range(len(file_names)):
    s = MolecularSystem()
    s.Read(MOL2File(),file_names[i][0])

    lt = MoltemplateLTFile()
    lt.AssignForcefieldType(s.molecules[0],file_names[i][1])

    # First, write the with-bond version
    lt.writeBondInfo = True
    with open("{}".format(file_names[i][2]),'w') as file:
        output.setoutput(file)
        s.Write(lt)

    # Then write the 'no bond' version
    lt.writeBondInfo = False
    with open("{}".format(file_names[i][3]),'w') as file:
        output.setoutput(file)
        s.Write(lt)
