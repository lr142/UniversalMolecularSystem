# Generate Special Structures such as Single Wall Carbon Nano Tube (SWCNT) or graphene sheet
from UniversalMolecularSystem import *
from BondDetection import *
import math
def SingleWallCarbonNanoTube(radius, length):
    # Returns a MolecularSystem containing a single molecule
    # the generated nanotube is parallel to the x-axis, with its left end centered at the origin. Translate and rotate
    # it if necessary.
    # Besides, this molecular system will be periodic.
    nanoTube = MolecularSystem()
    nanoTube.molecules.append(Molecule())
    nanoTube.molecules[0].name = "SWCNT"

    bondLength = 1.42 # The C-C bond length in graphene. Can be tuned to be more accurate
    A = bondLength / 2.0                    # A is the distance between site 1 and 2 in x-direction
    B = 2.0 * bondLength                    # B is the distance between site 1 and site 4 in a C6 ring
    C = bondLength * math.sqrt(3.0)         # C is the distance between parallel edges of a C6 ring


    # Radial Number of repeating units will be the closest integer. Therefore, each aromatic ring in the radial direction may be
    # slightly stretched or squeezed. The larger the tube, the smaller this effect is.
    NRadial = round(  math.pi / math.asin(C/2.0/radius)    )
    # The above strange equation is from some geometric calc:
    # Let theta = 2Pi/N  be the angle spanned by a ring
    # Then 2R*sin(theta/2) should equal the distance between 2 edges in a ring, ie. 2R*sin(Pi/N) = C

    # Number of repeating periodic units in the axial direction. The structure will be slightly larger than the given
    # length to make it periodic
    NAxial = int(math.ceil(length/(bondLength+B)))

    nanoTube.boundary = [[0,0,0],[0,0,0],[0,0,0]]
    nanoTube.boundary[0][0] = NAxial * (bondLength+B)
    nanoTube.boundary[1][1] = nanoTube.boundary[2][2] = 2 *(radius + 3.5)

    CAtom = Atom()
    CAtom.name = CAtom.type = 'C.ar'
    CAtom.element = 'C'
    CAtom.charge = 0.0
    for i in range(NAxial):
        for j in range(NRadial):
            # There are only four atoms in each primitive unit
            # specified in cylindrical coordination
            Xoffset = i * (bondLength + B)
            Xcoord = [Xoffset, Xoffset+A, Xoffset+A+bondLength, Xoffset+B ]
            Phioffset = j * (math.pi*2.0/NRadial)
            halfAngle = math.pi/NRadial
            Phicoord = [Phioffset, Phioffset+halfAngle, Phioffset+halfAngle, Phioffset]
            for n in range(4):
                atom = CAtom.Copy()
                atom.x = Xcoord[n]
                atom.y = radius * math.cos(Phicoord[n])
                atom.z = radius * math.sin(Phicoord[n])
                nanoTube.molecules[0].atoms.append(atom)

    # Find bonds
    nanoTube.RenumberAtomSerials()
    nanoTube.AutoDetectBonds(DefaultBondRules(2.0))

    return nanoTube


def Test():
    nt = SingleWallCarbonNanoTube(20,30)
    nt.Summary()
    with open('dump.mol2','w') as file:
        output.setoutput(file)
        nt.Write(MOL2File())
        output.setoutput(None)

if __name__ == '__main__':
    Test()



