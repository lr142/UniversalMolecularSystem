from Utility import *
from UniversalMolecularSystem import *
from MOL2File import *

# In every LAMMPS data file, the following info can be read:
# A comment line (including timestep) eg. LAMMPS data file via write_data, version 3 Mar 2020, timestep = 10000
# number of atoms, bonds, angle, dihedrals, etc.
# boundary info, e.g:
# -2.0000000000000000e+01 1.7000000000000000e+02 ylo yhi
# -2.0000000000000000e+01 1.7000000000000000e+02 zlo zhi
# Then the following sections:
# Masses
# Pair Coeffs # lj/cut/coul/long
# Bond Coeffs # harmonic
# Angle Coeffs # harmonic
# Dihedral Coeffs # opls
# Improper Coeffs # harmonic
# Atoms # full
# Velocities
# Bonds
# Angles
# Not all sections may be read. For example velocities may be absent.
# In this function we shall read in only boundary, masses, atoms, bonds.


# No element type is provided, but can be inferred from masses.
# The standard atomic weights are stored in 'atomicweights.csv' and are read only a
# LAMMPSDATAFile is created. Some atomic weights may not be included in the list,
# such as some united-atom. The user can call 'LAMMPSDATAFile.SetAtomicWeight(name, weight)'
# to set or override atomic weights.
# Matching criteria is withing [-0.1, +0.1]. types without a match shall be assigned element 'M'

class LAMMPSDATAFile(MolecularFile):
    element_weights: map
    atomtype_to_element: map
    globalSerialToLocalSerialMap: map
    globalSerialToMoleculeMap: map
    handleImageFlag: bool   # Whether to take the periodic image flag into consideration. Default is yes.

    def __init__(self,handleImageFlag = True):
        self.element_weights = {}
        import os
        default_route = '{}/{}'.format(os.path.dirname(__file__), "atomicweights.csv")
        # the first part of the path is needed if case the script is called from other directories
        with open(default_route,'r') as file:
            lines = file.readlines()
        for l in lines:
            parts = l.strip().split(',')
            if len(parts) == 2 and parts[1] != '':
                self.element_weights[parts[0]] = float(parts[1])

        self.atomtype_to_element = {}
        self.globalSerialToLocalSerialMap = {}
        self.globalSerialToMoleculeMap = {}
        self.handleImageFlag = handleImageFlag

    def SetAtomicWeight(selfs,element,weight):
        self.element_weights[element] = weight

    def __searchTypeFromWeight(self,weight):
        import math
        for t in self.element_weights:
            if math.fabs(weight - self.element_weights[t]) < 0.1:
                return t
        return 'M' # return M on non-found


    def Read(self,molecularSystem,filename):
        contents = None
        try:
            with open(filename,'r') as file:
                contents = file.readlines()
        except:
            error("Can't open LAMMPS DATA file {}".format(filename),False)
            return False

        def __Read_Masses(i):
            while(i < len(contents)):
                line = contents[i].strip()
                if line == '':
                    pass
                elif line[0].isalpha():
                    return i
                else:
                    parts = line.split()
                    self.atomtype_to_element[parts[0]] = self.__searchTypeFromWeight(float(parts[1]))
                i += 1

        def __Read_Atoms(i):
            # Temporary Atom List

            atoms_list = []
            maxMolID = 0
            while(i < len(contents)):
                line = contents[i].strip()
                if line == '':
                    pass
                elif line[0].isalpha():
                    break
                else:
                    parts = line.split()
                    atom = Atom()
                    atom.systemwideSerial = parts[0]
                    # Use the parent flag to store its molecule ID
                    atom.parent = parts[1]
                    maxMolID = max(maxMolID, int(atom.parent))
                    atom.element = atom.name = self.atomtype_to_element[parts[2]]
                    atom.type = "{}.{}".format(atom.element, parts[2])
                    atom.charge = float(parts[3])
                    atom.x = float(parts[4])
                    atom.y = float(parts[5])
                    atom.z = float(parts[6])
                    if self.handleImageFlag:
                        # Read periodic image flags. There may be flags in file. In which case there's nothing to do.
                        try:
                            imagex = int(parts[7])
                            imagey = int(parts[8])
                            imagez = int(parts[9])
                            atom.x += imagex * molecularSystem.boundary[0][0]
                            atom.y += imagey * molecularSystem.boundary[1][1]
                            atom.z += imagez * molecularSystem.boundary[2][2]
                        except:
                            pass

                    atoms_list.append(atom)

                i += 1
            # post processing:
            molecules = [ Molecule() for iMol in range(maxMolID)]
            for iMol, mol in enumerate(molecules,start=1):
                mol.serial = iMol
            for atom in atoms_list:
                iMol = int(atom.parent)
                belongToMol = molecules[iMol-1]
                atom.serial = '{}'.format(len(belongToMol.atoms)+1)
                self.globalSerialToLocalSerialMap[atom.systemwideSerial] = atom.serial
                self.globalSerialToMoleculeMap[atom.systemwideSerial] = belongToMol
                belongToMol.atoms.append(atom)

            molecularSystem.molecules = molecules
            return i

        def __Read_Bonds(i):
            while(i < len(contents)):
                line = contents[i].strip()
                if line == '':
                    pass
                elif line[0].isalpha():
                    return i
                else:
                    parts = line.split()
                    bond = Bond()
                    bond.type = parts[1]
                    b_from = parts[2]
                    b_to = parts[3]
                    mol1 = self.globalSerialToMoleculeMap[b_from]
                    mol2 = self.globalSerialToMoleculeMap[b_to]

                    if mol1.serial == mol2.serial:
                        bond.atom1 = self.globalSerialToLocalSerialMap[b_from]
                        bond.atom2 = self.globalSerialToLocalSerialMap[b_to]
                        mol1.bonds.append(bond)
                    else:
                        bond.atom1 = b_from
                        bond.atom2 = b_to
                        molecularSystem.interMolecularBonds.append(bond)
                i += 1

        def __Read_Boundary(i):
            # Assuming Orthogonal Box, and only read in size and origin of the box
            molecularSystem.boundary = [[0,0,0],[0,0,0],[0,0,0]]
            molecularSystem.origin = [0,0,0]
            for ixyz in range(3):
                parts = contents[i+ixyz].strip().split()
                lo = float(parts[0])
                hi = float(parts[1])
                molecularSystem.boundary[ixyz][ixyz] = hi-lo
                molecularSystem.origin[ixyz] = lo
            return i+3


        i = 0
        while(i < len(contents)):
            line = contents[i]
            if line.startswith('Masses'):
                i = __Read_Masses(i+1)
                # for t in self.atomtype_to_element:
                #     print('{} - {}'.format(t,self.atomtype_to_element[t]))
            elif line.startswith('Atoms'):
                i = __Read_Atoms(i+1)
            elif line.startswith('Bonds'):
                i = __Read_Bonds(i+1)
            elif i < 50 and line.find('hi') != -1:   # boundary info at the beginning
                i = __Read_Boundary(i)
            else:
                i = i+1


