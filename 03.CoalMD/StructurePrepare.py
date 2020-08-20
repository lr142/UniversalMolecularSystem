import sys
sys.path.append("../../")
from UniversalMolecularSystem import *
from Utility import *
from XSDFile import *
from PDBFile import *
from XYZFile import *
from BondDetection import *
from SpecialStructure import *
import os
import math



class StructruePrepParameter:
    # This is a a struct containing the following parameters that should be set by user before called 'StructurePrep'
    numberOfCoalCellsInXDirection: int
    extraPaddingInXDirection: float
    outerDiameter: float
    innerDiameter: float
    lengthOfFilledCH4: float   # fill the coal tube with CH4 in range [0, lengthOfFilledCH4)
    lengthOfFilledH2O: float
    # fill the coal tube with H2O in range [lengthOfFiledH [lengthOfFilledCH4, lengthOfFilledCH4 + lengthOfFilledH2O)

    # Following parameters are adjustable but given by default
    coalCellFile: str   # XSD file that comes also with cell size ( 99.7476 Å )
    methaneCellFile: str  # PDB file of CH4 in a cell ( 100 Å )
    waterCellFile: str # PDB file of H2O in a cell (100 Å)
    minDistanceForClashing: float # minimal distance (between heavy atoms) that is considered to be clashing

    def __init__(self):
        import os
        commonPath = '{}/{}'.format(os.path.dirname(__file__), "000.Common.Files.Read.Only")
        self.coalCellFile = '{}/{}'.format(commonPath, "CoalCell.xsd" )
        self.methaneCellFile = '{}/{}/{}'.format(commonPath, "water.and.ch4.cells", "4158methane.pdb")
        self.waterCellFile = '{}/{}/{}'.format(commonPath, "water.and.ch4.cells", "33428water.pdb")
        self.minDistanceForClashing = 3.0

def __CuttingTube(cell, length, outerDiameter, innerDiameter, centerY, centerZ):
    # Making a tube that is _length Å long (periodic), with outer shell diameter _outer Å, inner shell diameter _inner Å
    # The tube is parallel to the x-direction with its axis located at _centerY and _centerX.
    # Make sure that the to be cut coalCell is large enough! If the cell is not large enough, the cutting will have no effect
    # Requirement: NO INTERMOLECULAR BONDS! OR THERE MAY BE DANGLING BONDS AFTER CUTTING!
    #              Assuming the origin of the cell is at [0,0,0]
    # Note: Setting _inner diameter to be 0.0 can cut out a cylinder rather than a cube.
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

def _S1_PrepareCoal(para):
    primitiveCoalCell = MolecularSystem('CoalTube')
    primitiveCoalCell.Read(XSDFile(), para.coalCellFile)

    imagesX = para.numberOfCoalCellsInXDirection
    imagesY = imagesZ = math.ceil(para.outerDiameter / primitiveCoalCell.boundary[0][0])
    images = []
    for i in range(imagesX):
        for j in range(imagesY):
            for k in range(imagesZ):
                images.append([i, j, k])

    # Coordinates in this file are given in fractional coordinates, which is a little annoying
    primitiveCoalCell.FractionalToCartesianCoordinates()

    # The as read structure has some inter-molecular bonds (Some coal mols are splitted into halves, so can't just
    # delete those intermolecular bonds). Trying to eliminate them.
    primitiveCoalCell = ReduceSystemToOneMolecule(primitiveCoalCell)
    # This will preserve the 'boundary' info stored in primitiveCoalCell
    tempMS = BreakupMoleculeByConnectivity(primitiveCoalCell.molecules[0])
    primitiveCoalCell.molecules = tempMS.molecules

    largeCoalCell = primitiveCoalCell.Copy()
    DuplicateSystemPeriodically(largeCoalCell, images)

    # This value should be set in periodic system.
    largeCoalCell.boundary[0][0] *= imagesX


    __CuttingTube(largeCoalCell,
                    length=999999,  # We are not cutting the cube in the X-direction but leave it as is.
                    outerDiameter=para.outerDiameter,
                    innerDiameter=para.innerDiameter,
                    centerY=para.outerDiameter/2.0,
                    centerZ=para.outerDiameter/2.0)

    # If no padding, the cell is periodic in the x-direction. In this case, move all atoms with x-coords < 0 to
    # the right end of the cell
    if math.fabs(para.extraPaddingInXDirection) < 1E-3:
        A = largeCoalCell.boundary[0][0]
        for m in largeCoalCell.molecules:
            for a in m.atoms:
                if a.x <= 0:
                    a.x += A
                if a.x >= A:
                    a.x -= A


    return largeCoalCell

def __FillCoalWithSmallMol(coalCell,smallMolCell,XFrom,XTo,outerDiameter,centerY,centerZ,minDist,perMolecule = True):
    # SmallMol refers to Carbon Nanotube C atoms, CH4, H2O, etc.
    # When the perMolecule flag is set, the function only checks the coordinates of the 1st atom in the molecule.
    # When the 1st atom of the molecule doesn't clash with the coal, the whole molecule is added to the system.
    # This is the case if the smallmol is CH4 or H2O and the 1st atom is always the heavy atom ('H' or 'O')
    # When the perMolecule flag is False, the function checks coordinates of each atom. Atoms that are not clashing are
    # added into the system as a single molecule. In the later case the bonding info stored in the original mol will be
    # discarded.

    length = XTo - XFrom
    if math.fabs(length) < 1E-6:  # Return on zero length
        return

    imagesX = math.ceil(length / smallMolCell.boundary[0][0])
    imagesY = math.ceil( (outerDiameter/2.0 + centerY) / smallMolCell.boundary[1][1])
    imagesZ = math.ceil( (outerDiameter/2.0 + centerZ) / smallMolCell.boundary[2][2])


    images = []
    for i in range(imagesX):
        for j in range(imagesY):
            for k in range(imagesZ):
                images.append([i,j,k])

    DuplicateSystemPeriodically(smallMolCell,images)
    __CuttingTube(smallMolCell,
                length=length,
                outerDiameter=outerDiameter,  # inner diameter of coal =  outer diameter of filled CH4
                innerDiameter=0,
                centerY=centerY,
                centerZ=centerZ)

    smallMolCell.Translate(XFrom,0,0)

    # Now we need to remove those molecules that clash with coalCell. For this purpose we need to construct a
    # neighbor list for the coalCell
    atomsInCoal = []
    atomsToBeAdded = []
    for m in coalCell.molecules:
        atomsInCoal.extend(m.atoms)

    if perMolecule:
        # To save time, test only distances between heavy atoms.
        for m in smallMolCell.molecules:
            atomsToBeAdded.append(m.atoms[0])   # atoms[0] is either 'C' in CH4, or 'O' in H2O
    else:
        for m in smallMolCell.molecules:
            for a in m.atoms:
                atomsToBeAdded.append(a)


    grid = minDist * 1.1
    sys.stdout.write("Building a neighbor list with grid size {} for all atoms in the coal tube...".format(grid))
    coalNeighList = NeighborList(atomsInCoal,grid)
    sys.stdout.write('Done.\n')
    sys.stdout.write("Checking for clashes between the coalCell and atoms to be filled...")
    clashing = coalNeighList.DetectClashing(atomsToBeAdded,minDist)


    removed = 0

    if perMolecule:
        newMolList = []
        for i,m in enumerate(smallMolCell.molecules):
            if not clashing[i]:
                newMolList.append(m)
            else:
                removed += 1
        sys.stdout.write("Done. {} mols excluded due to clashing, {} mols added\n".format(removed,len(newMolList)))
        smallMolCell.molecules = newMolList
    else:
        newMol = Molecule()
        for i,a in enumerate(atomsToBeAdded):
            if not clashing[i]:
                newMol.atoms.append(a)
            else:
                removed += 1
        sys.stdout.write("Done. {} atoms excluded due to clashing, {} atoms added\n".format(removed,len(newMol.atoms)))
        smallMolCell.molecules = [newMol]

    ExtendSystem(coalCell,smallMolCell)


def _S2_1_PaddingWithCarbonNanoTubeWall(coalCell,para):
    # Pad the inner wall of the coal tube with a double-walled carbon nanotube
    nanoTube = SingleWallCarbonNanoTube(para.innerDiameter/2.0 + 3 * para.minDistanceForClashing, coalCell.boundary[0][0])
    wall2 = SingleWallCarbonNanoTube(para.innerDiameter/2.0 + para.minDistanceForClashing, coalCell.boundary[0][0])

    nanoTube.molecules[0].atoms.extend(wall2.molecules[0].atoms)
    nanoTube.RenumberAtomSerials()

    # Alignment
    nanoTube.Translate(0, para.outerDiameter / 2, para.outerDiameter / 2)

    # In order to check for clashing, split the nanoTube to individual atoms, i.e. 1 atom per molecule
    newMols = []
    for a in nanoTube.molecules[0].atoms:
        mol = Molecule()
        mol.name = nanoTube.molecules[0].name
        mol.atoms.append(a.Copy())
        newMols.append(mol)

    nanoTube.molecules = newMols
    nanoTube.RenumberAtomSerials()
    __FillCoalWithSmallMol(coalCell,nanoTube,
                           XFrom = 0.0,
                           XTo = coalCell.boundary[0][0] - para.minDistanceForClashing * 0.5,
                           outerDiameter = para.innerDiameter + 6 * para.minDistanceForClashing,
                           centerY = para.outerDiameter/2.0,
                           centerZ = para.outerDiameter/2.0,
                           minDist = para.minDistanceForClashing * 0.5,  # We want a 'tight' wall that doesn't leak
                           perMolecule = False
                           )


def _S2_2_FillCoalWithMethane(coalCell,para):

    if math.fabs(para.lengthOfFilledCH4) < 1E-6:
        sys.stdout.write('Skipped.\n')
        return

    methaneCell = MolecularSystem('MethaneMols')
    methaneCell.Read(PDBFile(),para.methaneCellFile)
    methaneCell.boundary = [[100,0,0],[0,100,0],[0,0,100]]
    methaneCell.AutoDetectBonds(DefaultBondRules(1.2))

    __FillCoalWithSmallMol(coalCell,methaneCell,
                           XFrom= 0.0,
                           XTo = para.lengthOfFilledCH4,
                           outerDiameter = para.innerDiameter,  # The outer D of filled mols is the inner D of the tube
                           centerY = para.outerDiameter/2.0,
                           centerZ = para.outerDiameter/2.0,
                           minDist = para.minDistanceForClashing)



def _S2_3_FillCoalWithWater(coalCell,para):

    if math.fabs(para.lengthOfFilledH2O) < 1E-6:
        sys.stdout.write('Skipped.\n')
        return

    waterCell = MolecularSystem('WaterMols')
    waterCell.Read(PDBFile(), para.waterCellFile)
    waterCell.boundary = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
    waterCell.AutoDetectBonds(DefaultBondRules(1.2))

    __FillCoalWithSmallMol(coalCell, waterCell,
                           XFrom= para.lengthOfFilledCH4,
                           XTo = para.lengthOfFilledCH4 + para.lengthOfFilledH2O,
                           outerDiameter=para.innerDiameter,  # The outer D of filled mols is the inner D of the tube
                           centerY=para.outerDiameter / 2.0,
                           centerZ=para.outerDiameter / 2.0,
                           minDist = para.minDistanceForClashing)

def __WriteTheSWCNTWallLTFile(wallMol,filename):
    # The format is:
    #import "oplsaa.lt"
    # Given1 inherits OPLSAA{
    # write('Data Atoms'){
    #      $atom:C1   $mol:.   @atom:92   0.0    -8.96237   -3.16977   -0.319371
    #      ......
    # }
    lines = []
    lines.append('import "oplsaa.lt"')
    lines.append('SWCNTWall inherits OPLSAA{')
    lines.append("  write('Data Atoms'){")
    for i,a in enumerate(wallMol.atoms,start=1):
        lines.append("     $atom:C{}   $mol:.   @atom:92   0.0    {}   {}   {}".format(i,a.x,a.y,a.z))
    lines.append('  }')
    lines.append('}')

    with open(filename,'w') as file:
        output.setoutput(file)
        for l in lines:
            output(l)
        output.setoutput(None)

def _S3_SortByMolAndGenerateSystemLT(cell,para):

    g1mols = []
    g2mols = []
    fsmols = []
    ch4mols = []
    watmols = []
    wallmol = None

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
            m.type = "Given2"
            g2mols.append(m)
            coalAtomCount += NAtoms
        elif NAtoms == 250:
            m.type = "FuchsSandoff"
            fsmols.append(m)
            coalAtomCount += NAtoms
        elif NAtoms == 5:
            m.type = "CH4"
            ch4mols.append(m)
            ch4AtomCount += NAtoms
        elif NAtoms == 4:
            m.type = "WAT"
            watmols.append(m)
            watAtomCount += NAtoms
        else:
            # Check if it's the Wall molecule:
            ifWallMol = True
            for a in m.atoms:
                if a.element.upper() != 'C':
                    ifWallMol = False
                    break

            if ifWallMol:
                m.type = "WALL"
                wallmol = m
                coalAtomCount += len(wallmol.atoms)

            else:
                error("Molecule with {} atoms should not exist in this system. Double check the preparation process.".format(NAtoms))

    newMolList = []
    newMolList.extend(g1mols)
    newMolList.extend(g2mols)
    newMolList.extend(fsmols)
    if wallmol != None:
        newMolList.append(wallmol)
    newMolList.extend(ch4mols)
    newMolList.extend(watmols)

    cell.molecules = newMolList
    cell.RenumberAtomSerials()

    # The output sequence shall be:
    # 1. The CoalTube itself
    # 2. SWCNT wall. This requires a separate custom-made LT file.
    # 3. Methane
    # 4. Water
    # 1 and 2 forms group gCoal
    # 3 and 4 constituent groups gCH4 and gH2O, respectively.
    # 3 and 4 together are group gMobile

    #output.setoutput()

    lines = [
        'import "../common.header.lt"',
        'write_once("Data Boundary"){',
        ]
    xmin = 0.0
    xmax = cell.boundary[0][0]*para.numberOfCoalCellsInXDirection + para.extraPaddingInXDirection
    ymin = zmin = -20.0
    ymax = zmax = para.outerDiameter + 20.0

    lines.append('     {} {} xlo xhi'.format(xmin,xmax))
    lines.append('     {} {} ylo yhi'.format(ymin,ymax))
    lines.append('     {} {} zlo zhi'.format(zmin,zmax))
    lines.append('}')
    lines.append('write_once("In Init"){')
    lines.append('     kspace_style pppm 0.001')
    lines.append('     boundary p p p')
    lines.append('}')
    lines.append('g1 = new Given1[{}]'.format(len(g1mols)))
    lines.append('g2 = new Given2[{}]'.format(len(g2mols)))
    lines.append('fs = new FuchsSandoff[{}]'.format(len(fsmols)))

    if wallmol != None:
        __WriteTheSWCNTWallLTFile(wallmol,"swcntwall.lt")
        lines.append('import "swcntwall.lt"')
        lines.append('wall = new SWCNTWall')

    if len(ch4mols) > 0:
        lines.append('ch4 = new Methane[{}]'.format(len(ch4mols)))
    if len(watmols) > 0:
        lines.append('h2o = new TIP4P[{}]'.format(len(watmols)))


    lines.append('write_once("In Settings"){')
    lines.append('     group gCoal id {}:{}'.format(1,coalAtomCount))
    if ch4AtomCount > 0:
        lines.append('     group gCH4 id {}:{}'.format(coalAtomCount+1, coalAtomCount+ch4AtomCount))
    if watAtomCount > 0:
        lines.append('     group gH2O id {}:{}'.format(coalAtomCount+ch4AtomCount+1, coalAtomCount+ch4AtomCount+watAtomCount))
    lines.append('     group gMobile id {}:{}'.format(coalAtomCount+1,coalAtomCount+ch4AtomCount+watAtomCount))
    lines.append('}')

    with open('system.lt','w') as file:
        output.setoutput(file)
        for l in lines:
            output(l)
        output.setoutput(None)


def Prepare(parameters):
    output('Step 1: Prepare the coal tube...')
    coalCell = _S1_PrepareCoal(parameters)
    coalCell.Summary()

    output('Step 2-1: Padding with a SWCNT...')
    _S2_1_PaddingWithCarbonNanoTubeWall(coalCell,parameters)
    coalCell.Summary()

    output('Step 2-2: Fill with CH4...')
    _S2_2_FillCoalWithMethane(coalCell,parameters)
    coalCell.Summary()

    output('Step 2-3: Fill with H2O...')
    _S2_3_FillCoalWithWater(coalCell,parameters)
    coalCell.Summary()

    output('Step 3: Finalize and generating files...')
    _S3_SortByMolAndGenerateSystemLT(coalCell,parameters)


    reducedCell = ReduceSystemToOneMolecule(coalCell)


    # for i in range(len(reducedCell.molecules[0].atoms)):
    #     if i%100 == 0:
    #         print(i)
    #     if reducedCell.molecules[0].atoms[i].element.upper() == 'H' or reducedCell.molecules[0].atoms[i].element.upper() == 'M':
    #         continue
    #     for j in range(i+1,len(reducedCell.molecules[0].atoms)):
    #         d = Distance(reducedCell.molecules[0].atoms[i],reducedCell.molecules[0].atoms[j])
    #         if d < 0.8:
    #             print('Collision between {} and {}'.format(i,j))


    reducedCell.Summary()
    with open('positions.xyz','w') as file:
        output.setoutput(file)
        reducedCell.Write(XYZFile())
        output.setoutput(None)




