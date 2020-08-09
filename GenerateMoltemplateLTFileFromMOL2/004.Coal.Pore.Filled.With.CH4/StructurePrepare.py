import sys
sys.path.append("../../")
from UniversalMolecularSystem import *
from Utility import *
from XSDFile import *
from PDBFile import *
from XYZFile import *


# Making a tube that is _length Å long (periodic), with outer shell diameter _outer Å, inner shell diameter _inner Å
# Make sure that the to be cut coalCell is large enough! If the cell is not large enough, the cutting will have no effect
# Requirement: NO INTERMOLECULAR BONDS! OR THERE MAY BE DANGLING BONDS AFTER CUTTING!
#              Assuming the origin of the cell is at [0,0,0]
def CuttingTube(coalCell,_length,_outer,_inner):
    molecule_remain = []
    def MolCenter(mol):
        # Geometric Center of molecules, not center of gravity
        cx = cy = cz = 0.0
        for a in mol.atoms:
            cx += a.x
            cy += a.y
            cz += a.z
        cx/=len(mol.atoms)
        cy/=len(mol.atoms)
        cz/=len(mol.atoms)
        return [cx,cy,cz]

    for mol in coalCell.molecules:
        center = MolCenter(mol)
        from math import pow
        from math import sqrt
        dist_to_axis = sqrt(pow(center[1] - _outer/2,2) + pow(center[2] - _outer/2,2))
        if center[0] > _length or dist_to_axis > _outer/2 or dist_to_axis < _inner/2:
            molecule_remain.append(False)
        else:
            molecule_remain.append(True)

    newMolList = []
    for i in range(len(molecule_remain)):
        if molecule_remain[i]:
            newMolList.append(coalCell.molecules[i])

    coalCell.molecules = newMolList
    coalCell.RenumberAtomSerials
    return True



def Prepare():
    images = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                images.append([i,j,k])

    primitiveCoalCell = MolecularSystem()
    primitiveCoalCell.Read(XSDFile(), 'CoalCell.xsd')
    # Coordinates in this file are given in fractional coordinates, which is a little annoying
    primitiveCoalCell.FractionalToCartesianCoordinates()

    primitiveCoalCell = ReduceSystemToOneMolecule(primitiveCoalCell)  # The as read structure has some inter-molecular bonds. Trying to elniminate them.
    tempMS = BreakupMoleculeByConnectivity(primitiveCoalCell.molecules[0])
    primitiveCoalCell.molecules = tempMS.molecules

    largeCoalCell = primitiveCoalCell.Copy()
    largeCoalCell.Summary()
    DuplicateSystemPeriodically(largeCoalCell,images)
    largeCoalCell.Summary()


    CuttingTube(largeCoalCell,200,100,50)
    largeCoalCell.Summary()


    output.setoutput(open('dump.mol2','w'))
    largeCoalCell.Write(MOL2File())


Prepare()



