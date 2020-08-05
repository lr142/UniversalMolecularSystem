import sys

sys.path.append("..")

from Utility import *
from UniversalMolecularSystem import *
from MOL2File import *

class MoltemplateLTFile(MolecularFile):
    def __init__(self):
        pass

    def AssignForcefieldType(self,molecule,force_field_type_file):
        with open(force_field_type_file,"r") as file:
            lines = file.readlines()

        serialToAtomMap = {}
        for atom in molecule.atoms:
            serialToAtomMap[atom.serial] = atom   # saves links to corresponding atoms

        for l in lines:
            pars = l.strip().split(',')
            serial = pars[0]
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




# if len(sys.argv) == 1:
#     sys.argv = ["","Given1.csv"]

file_names = [["Given1.mol2","Given1.csv","given1.lt"],
            ["Given2.mol2","Given2.csv","given2.lt"],
            ["FuchsSandoff.mol2","FuchsSandoff.csv","fuchssandoff.lt"]
         ]

for i in range(3):
    s = MolecularSystem()
    s.Read(MOL2File(),file_names[i][0])
    s.molecules[0].Check(True) # Force renaming to remove duplicate atom names

    lt = MoltemplateLTFile()
    lt.AssignForcefieldType(s.molecules[0],file_names[i][1])

    file = open("moltemplate/{}".format(file_names[i][2]),'w')
    output.setoutput(file)
    # s.Write(MOL2File())

    s.Write(lt)



