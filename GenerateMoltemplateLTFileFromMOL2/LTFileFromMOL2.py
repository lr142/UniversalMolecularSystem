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

        for l in lines:
            output(l)


if len(sys.argv) == 1:
    sys.argv = ["","Given1.csv"]

lt = MoltemplateLTFile()
lt.AssignForcefieldType(None,sys.argv[1])


