# Adjust for TIP4P water. Search for water molecules in the system; Insert an 'M' atom between O and 2H atoms, and
# adjust the O-H bond length and H-O-H bond angle so that it equals values specified in TIP4P forcefield:
# O-H length = 0.95721 Å, H-O-H angle = 104.5183 degrees.
#         $atom:O      $mol:.    @atom:65     0.0      0.0    0.0     0.0
#         $atom:M      $mol:.    @atom:67     0.0      0.0    0.0     0.150
#         $atom:H1     $mol:.    @atom:66     0.0      0.0    -0.75695     0.5859
#         $atom:H2     $mol:.    @atom:66     0.0      0.0    +0.75695     0.5859

from MolecularManipulation import *
from UniversalMolecularSystem import *
from Utility import *

def PatchTIP4P(mol):

    # Let us define an axis system like this:
    # -----------O------------->   y
    #          / | \
    #         /  M  \
    #       H1   |   H2
    #            v
    #            z
    # x axis is perpendicular to the viewer and pointing inwards.


    aM = Atom()
    aM.element = aM.name = "M"
    mol.atoms.insert(1,aM)

    aO = mol.atoms[0]
    aH1 = mol.atoms[2]
    aH2 = mol.atoms[3]

    vO = [aO.x,aO.y,aO.z]

    vZ = VectorMinus([(aH1.x+aH2.x)/2, (aH1.y+aH2.y)/2, (aH1.z+aH2.z)/2],vO)
    vZ = VectorNormalize(vZ)

    vOH1 = VectorMinus([aH1.x,aH1.y,aH1.z], vO)
    vX = VectorCrossProduct(vZ, vOH1)
    vX = VectorNormalize(vX)

    vY = VectorCrossProduct(vZ, vX)

    # Calculate the positions of atoms in the new system.
    # vMid is the position of the mid point between 2 Hs relative to O
    vOM = VectorScalarProduct(vZ, 0.150)
    vMid = VectorScalarProduct(vZ, 0.5859)
    vOH1 = VectorAdd(vMid, VectorScalarProduct(vY,-0.75695))
    vOH2 = VectorAdd(vMid, VectorScalarProduct(vY,+0.75695))

    vM = VectorAdd(vO, vOM)
    vH1 = VectorAdd(vO, vOH1)
    vH2 = VectorAdd(vO, vOH2)

    aM.x, aM.y, aM.z = vM
    aH1.x, aH1.y, aH1.z  = vH1
    aH2.x, aH2.y, aH2.z = vH2



def AdjustForTIP4P(ms):
    from BondDetection import DefaultBondRules
    ms.AutoDetectBonds(DefaultBondRules(2.5),)
    newMS = BreakupMoleculeByConnectivity(ms.molecules[0])
    patched = False

    for m in newMS.molecules:
        if len(m.atoms) == 3 \
            and m.atoms[0].element.upper() == 'O' \
                and m.atoms[1].element.upper() == 'H' \
                and m.atoms[2].element.upper() == 'H':
            PatchTIP4P(m)
            patched = True

    if patched:
        newMS = ReduceSystemToOneMolecule(newMS)
        ms.molecules = newMS.molecules
        ms.RenumberAtomSerials()





