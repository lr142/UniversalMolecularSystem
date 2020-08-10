import sys
sys.path.append("../../")
from UniversalMolecularSystem import *
from Utility import *
from XSDFile import *
from PDBFile import *
from XYZFile import *
from BondDetection import *
import os


# Making a tube that is _length Å long (periodic), with outer shell diameter _outer Å, inner shell diameter _inner Å
# The tube is parallel to the x-direction with its axis located at _centerY and _centerX.
# Make sure that the to be cut coalCell is large enough! If the cell is not large enough, the cutting will have no effect
# Requirement: NO INTERMOLECULAR BONDS! OR THERE MAY BE DANGLING BONDS AFTER CUTTING!
#              Assuming the origin of the cell is at [0,0,0]
# Note: Setting _inner diameter to be 0.0 can cut out a cylinder rather than a cube.
def CuttingTube(cell, length, outerDiameter, innerDiameter, centerY, centerZ):
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

    for mol in cell.molecules:
        center = MolCenter(mol)
        from math import pow
        from math import sqrt
        dist_to_axis = sqrt(pow(center[1] - centerY, 2) + pow(center[2] - centerZ, 2))
        if center[0] > length or dist_to_axis > outerDiameter/2 or dist_to_axis < innerDiameter/2:
            molecule_remain.append(False)
        else:
            molecule_remain.append(True)

    newMolList = []
    for i in range(len(molecule_remain)):
        if molecule_remain[i]:
            newMolList.append(cell.molecules[i])

    cell.molecules = newMolList
    cell.RenumberAtomSerials()
    return True

def PrepareCoal():
    images = []
    for i in range(5):
        for j in range(2):
            for k in range(2):
                images.append([i,j,k])

    primitiveCoalCell = MolecularSystem('CoalTube')
    primitiveCoalCell.Read(XSDFile(), 'CoalCell.xsd')
    # Coordinates in this file are given in fractional coordinates, which is a little annoying
    primitiveCoalCell.FractionalToCartesianCoordinates()

    # The as read structure has some inter-molecular bonds (Some coal mols are splitted into halves, so can't just
    # delete those intermolecular bonds). Trying to eliminate them.
    primitiveCoalCell = ReduceSystemToOneMolecule(primitiveCoalCell)
    # This will preserve the 'boundary' info stored in primitiveCoalCell
    tempMS = BreakupMoleculeByConnectivity(primitiveCoalCell.molecules[0])
    primitiveCoalCell.molecules = tempMS.molecules

    largeCoalCell = primitiveCoalCell.Copy()
    largeCoalCell.Summary()
    DuplicateSystemPeriodically(largeCoalCell,images)
    largeCoalCell.Summary()


    CuttingTube(largeCoalCell, length=500, outerDiameter=100, innerDiameter=50, centerY=50, centerZ=50)
    largeCoalCell.Summary()

    return largeCoalCell

def FillCoalWithMethane(coalCell):
    methaneCell = MolecularSystem('MethaneMols')
    methaneCell.Read(PDBFile(),"4158methane.pdb")
    methaneCell.boundary = [[100,0,0],[0,100,0],[0,0,100]]
    methaneCell.AutoDetectBonds(DefaultBondRules(1.2))

    images = []
    for i in range(4):
        images.append([i,0,0])
    DuplicateSystemPeriodically(methaneCell,images)
    CuttingTube(methaneCell,length=400,outerDiameter=50,innerDiameter=0,centerY=50,centerZ=50)

    # Now we need to remove those CH4 molecules that clash with coalCell. For this purpose we need to construct a
    # neighbor list for the coalCell
    atomsInCoal = []
    atomsCinCH4 = []
    for m in coalCell.molecules:
        atomsInCoal.extend(m.atoms)
    # To save time, test only distances between C atoms in CH4 and coal atoms. For this purpose, we need a list of all C atoms in CH4.
    for m in methaneCell.molecules:
        atomsCinCH4.append(m.atoms[0])   # atoms[0] is always C.

    methaneCell.Summary()
    output("Building a neighbor list for all atoms in the coal tube...")
    coalNeighList = NeighborList(atomsInCoal,3.0)
    output("Done.")
    output("Checking for clashing between CH4 and the coalCell...")
    clashing = coalNeighList.DetectClashing(atomsCinCH4,3.0)
    output("Done.")
    newMolList = []
    for i,m in enumerate(methaneCell.molecules):
        if clashing[i] == False:
            newMolList.append(m)
    methaneCell.molecules = newMolList
    methaneCell.Summary()

    ExtendSystem(coalCell,methaneCell)

def SortByMolAndGenerateSystemLT(cell,path):

    g1mols = []
    g2mols = []
    fsmols = []
    ch4mols = []
    watmols = []

    coalAtomCount = 0
    ch4AtomCount = 0
    watAtomCount = 0

    for m in cell.molecules:
        NAtoms = len(m.atoms)
        if NAtoms == 130:
            m.type = "Given1"
            g1mols.append(m)
            coalAtomCount += NAtoms
        elif NAtoms == 188:
            m.type == "Given2"
            g2mols.append(m)
            coalAtomCount += NAtoms
        elif NAtoms == 250:
            m.type == "FuchsSandoff"
            fsmols.append(m)
            coalAtomCount += NAtoms
        elif NAtoms == 5:
            m.type == "CH4"
            ch4mols.append(m)
            ch4AtomCount += NAtoms
        elif NAtoms == 3:
            m.type == "WAT"
            watmols.append(m)
            watAtomCount += NAtoms
        else:
            error("Molecule with {} atoms should not exist in this system. Double check the preparation process.".format(NAtoms))

    newMolList = []
    newMolList.extend(g1mols)
    newMolList.extend(g2mols)
    newMolList.extend(fsmols)
    newMolList.extend(ch4mols)
    newMolList.extend(watmols)

    cell.molecules = newMolList
    cell.RenumberAtomSerials()

    output.setoutput(open('{}/system.lt'.format(path),'w'))
    output('import "oplsaa.lt"\n'
           'import "../given1.nb.lt"\n'
           'import "../given2.nb.lt"\n'
           'import "../fuchssandoff.nb.lt"\n'
           'import "../methane.lt"\n'
           'import "../tip4p.lt"\n'
           'import "../tip4p2005.lt"\n'
           'write_once("Data Boundary"){\n')
    output('     {} {} xlo xhi'.format(0,cell.boundary[0][0]*5))
    output('     -15 115 ylo yhi\n'
           '     -15 115 zlo zhi\n'
           '}\n'
           'write_once("In Init"){\n'
           '     boundary p p p\n'
           '}')
    output('g1 = new Given1[{}]'.format(len(g1mols)))
    output('g2 = new Given2[{}]'.format(len(g2mols)))
    output('fs = new FuchsSandoff[{}]'.format(len(fsmols)))
    if len(ch4mols) > 0:
        output('ch4 = new Methane[{}]'.format(len(ch4mols)))
    if len(watmols) > 0:
        output('ch4 = new TIP4P2005[{}]'.format(len(watmols)))

    if len(ch4mols) > 0:
        output('write_once("In Settings"){{\n'
               '     group gCH4 id {}:{}\n'
               '}}'.format(coalAtomCount+1, coalAtomCount+ch4AtomCount)
               )

    if len(watmols) > 0:
        output('write_once("In Settings"){{\n'
            '     group gWAT id {}:{}\n'
            '}}'.format(coalAtomCount+ch4AtomCount+1, coalAtomCount+ch4AtomCount+watAtomCount))


    output.setoutput(None)


def Main():
    largeCoalCell = PrepareCoal()

    FillCoalWithMethane(largeCoalCell)
    largeCoalCell.Summary()

    path = 'moltemplate.working'
    SortByMolAndGenerateSystemLT(largeCoalCell,path)

    largeCoalCell = ReduceSystemToOneMolecule(largeCoalCell)

    largeCoalCell.Summary()

    output.setoutput(open('{}/positions.xyz'.format(path),'w'))
    largeCoalCell.Write(XYZFile())


if __name__ == '__main__':
    Main()



