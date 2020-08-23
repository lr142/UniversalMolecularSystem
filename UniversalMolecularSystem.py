#!/usr/bin/python
import sys
import math
import numpy as np
from Utility import *


class Bond:
    atom1: str
    atom2: str
    type: str
    length: float
    def __init__(self):
        self.atom1 = None
        self.atom2 = None
        self.type = None
        self.length = 0.0

    def Copy(self):
        import copy
        return copy.copy(self)

class Atom:
    element: str
    x:float
    y:float
    z:float
    type:str
    charge:float
    flexible:bool
    serial:str
    systemwideSerial: str
    layerInfo: str
    parent: str

    # Note: If complex data structures just as lists or references to other objects are included in Atom in future
    # versions, remember that Atom.Copy() is used to duplicate an atom, which uses copy.copy() internally to
    # 'shallow' copy atoms. Make sure such operations are safe
    # and reasonable. For example, if an atom contains a link to a third object, keep in mind that the external object
    # will not be copied. This note applies also to 'Bond'.
    # Since both Atom and Bond are considered to be 'elementary' classes, try not to include references or complex stuffs in them.

    def __init__(self):
        self.element = None   # Its element, must be correctly capitalized like "C" or "Rh"; "C1", "c", "rh" are not acceptable.
        self.x = self.y = self.z = 0.0  # Cartesian or fractional coordinates
        self.name = None   # Its internal name as appears in mol2 files, like "H13" or "C2". No need to be unique.
        self.type = None   # Its force field type or its hybridization type, like "C.3" or "P.2"
        self.charge = 0.0
        self.flexible = None
        self.serial = None     # serial number in a molecular, usually a number starting from 1 but can be any string
                               # serial number must be unique within a molecule
        self.systemwideSerial = None  # unique serial number in the molecular system. Must be present if there are
                               # intermolecular bonds (including hydrogen bonds) within the system.
        self.layerInfo = None  # Used in QM/MM models, ie. high/mid/low layer in ONIOM methods
        self.parent = None     # Used in some cases to record its parent, for example the molecule or residue it belongs to

        # Not all properties are available or relevant in all cases. File readers/writers and users can add additional
        # properties to an Atom if necessary.

    def Copy(self):
        # Returns a copy of the current Atom. It uses 'shallow' copy. See the notes at the beginning of this class
        import copy
        return copy.copy(self)

    def Translate(self,dx,dy,dz):  # Self-evident
        self.x += dx
        self.y += dy
        self.z += dz

    def ShowAsXYZ(self):
        print(str(self.element) + " " +
              str(self.x) + " " +
              str(self.y) + " " +
              str(self.z))
    def ShowAllFields(self):
        print(str(self.serial) + " " +
              str(self.element) + " " +
              str(self.x) + " " +
              str(self.y) + " " +
              str(self.z) + " " +
              str(self.name) + " " +
              str(self.type) + " " +
              str(self.flexible) + " " +
              str(self.charge)
              )

class Molecule:
    atoms: [Atom]
    bonds: [Bond]
    name: str
    bondedTo: [map]
    serial: str
    type: str

    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.name = None  # Names are user-defined identifiers, e.g. it can be type+serial in a protein. No need to be unique
        self.serial = None # Such as its serial in a peptide chain, must be unique within a MolecularSystem
        self.type = None # Such as the 3-Letter code residue name GLY
    def FindBonds(self,rules):
        for i in range(0,len(self.atoms)):
            #sys.stderr.write("Atom No. {}\n".format(i))
            for j in range(i+1,len(self.atoms)):
                a1 = self.atoms[i]
                a2 = self.atoms[j]
                length = math.sqrt (  (a1.x-a2.x)**2 +
                                 (a1.y-a2.y)**2 +
                                 (a1.z-a2.z)**2)
                order = rules.CheckRules(a1.name,a2.name,length)
                if ( order != None ):
                    bond = Bond()
                    bond.atom1 = int(i)
                    bond.atom2 = int(j)
                    bond.length = length
                    bond.type = order
                    self.bonds.append(bond)

    def BondedMap(self):
        # An auxiliary function, returns a list of sets (called bondedTo), for example
        # bondedTo[0] records the indexes (not serials) of atoms that are connected to atom with index (not serial) 0
        # Returns None if there are dangling bonds.
        bondedMap = [set() for i in range(len(self.atoms))]

        serialToIndexMap = {}
        for i,atom in enumerate(self.atoms):
            serialToIndexMap[atom.serial] = i

        for b in self.bonds:
            fromSerial = b.atom1
            toSerial = b.atom2
            if fromSerial in serialToIndexMap and toSerial in serialToIndexMap:
                fromIndex = serialToIndexMap[fromSerial]
                toIndex = serialToIndexMap[toSerial]
                bondedMap[fromIndex].add(toIndex)
                bondedMap[toIndex].add(fromIndex)
            else:
                return None

        return bondedMap

    def CheckConsistency(self):
        # If all checks are passed, return True, otherwise False
        # This is a consistency check to make sure that:
        # 1. Atom serial numbers must be unique
        serialToIndexMap = {}
        for i,atom in enumerate(self.atoms):
            serial = atom.serial
            if serial in serialToIndexMap:  # serial are not unique
                error("Atom serials in the molecule are not unique, found that serial {} appears again".format(serial),False)
                return False

            serialToIndexMap[serial] = i

        # 2. No dangling bonds, i.e no bonds pointing to atoms that are not exist.
        #    Used the auxiliary function 'BondedMap'
        result = self.BondedMap()
        if result == None:
            error("Dangling bond found in molecule: {} - {} with type {}".format(b.atom1,b.atom2,b.type),False)
            self.bondedTo = None  # Destroys the bondedMap and quit
            return False

        return True

    def Copy(self):
        import copy
        newMol = copy.copy(self)   # A 'Shallow' copy that copies only references to its string and integer properties
        newMol.atoms = []
        newMol.bonds = []
        for a in self.atoms:
            newMol.atoms.append(a.Copy())
        for b in self.bonds:
            newMol.bonds.append(b.Copy())
        return newMol


    def Summary(self):
        output("Molecule Serial: {}, Name: {}, Type: {}, with {} atoms, and {} bonds".format(
            self.serial,self.name,self.type,len(self.atoms),len(self.bonds)
        ))

class MolecularFile:
    # This is an abstract class that represents a molecular system description file and the methods to Read/Write
    # the file. For each supported file type, such as .xyz, .mol2, .pdb, etc., a separate concrete class shall be
    # created which implements the Read() and Write() virtual functions.
    def __init__(self):  # Virtual function
        pass
    def Read(self,molecularSystem,filename):   # Virtual function
        # Remember that any "Read" operation (even if not successful) will reset the molecular system
        # In cases where there is a need to read a system from multiple files, read them into separate
        # MolecularSystems and merge them.
        pass
    def Write(self,molecularSystem):   # Virtual function
        pass
class BondDetector:
    # This is an abstract class that has a sole function called Detect, which recognize and generates bonds within
    # a molecular system. In 'BondDetection.py' we implemented a clumsy detector based solely on atom distances and element
    # types. In the future we can write fancier detectors like one based on machine learning.
    def Detect(self,molecularSystem,flushCurrentBonds):
        pass

class Trajectory:
    # A trajectory is a component of a molecular system, generated by geometric optimization or MD or MC
    # The Trajectory class utilizes Numpy 'narray's to store positions, velocities, forces, etc of each atom in each frame
    timestep: float # timestep length in the unit of fs, (consistent with the 'real' units in LAMMPS)
    NFrames: int     # number of frames
    timesteps_of_each_frame: []  # records the timestep each frame corresponds to.
    NAtoms: int     # number of total atoms
    serial_to_index_map: map   # A map that maps atom.systemwideSerial to indexes in the following array
    index_to_serial: []        # with length NAtom, this array records each atom's systemwideSerial to its MolecularSystem
    # Below are several lists of narrays that records the coordinates, velocities, and forces of each atom in each frame.
    # Not all fields may be present in the dump file,for example velocities and forces may be absent (set to None)
    # These lists of narrays are in fact three-dimensional, with the highest (1st) dimension being the index of the frame,
    #  the 2nd dim being the atom index, and the last (3rd) dim being x,y,and z. Therefore, positions[100][25][2]
    # refers to the z coordinate of the 26th atom in the 101st frame.
    # NOTE: I tried to implement them as 3-dimensional narrays. But since we don't know the number of frame beforehand,
    # dynamically expanding the narrays is a great pain.
    positions: []
    velocities: []
    forces: []
    dtype: str    # data type for float numbers. Default 'float32' (double). Can specify 'float16' at creation to save space
    # parentMolecularSystem:   # Reference to its molecular system. May be None.

    # Assumptions:
    # All frames have the same number of atom. i.e: no lost atoms. The program does have the capability
    # to handle lost atoms, in which case the missing atom will be assigned all-zero xyz/velocity/force.

    def __init__(self, parentMolecularSystem = None, timestep_in_fs = 1.0, dtype = 'float32'):
        # it's the caller's responsibility to tell the program how long is each timestep, which can't be read from
        # LAMMPS dump files. A default value (1.0 fs) is given since in some cases ( like goemetry optimization in QM)
        # a timestep is irrelevant.
        self.timestep = timestep_in_fs
        self.timesteps_of_each_frame = []
        self.NFrames = 0
        self.NAtoms = 0
        self.dtype = dtype
        self.positions = self.velocities = self.forces = 0
        self.parentMolecularSystem = parentMolecularSystem
        self.serial_to_index_map = {}  # A map that maps atom.systemwideSerial to indexes in the following array
        self.index_to_serial = []
        self.positions = []
        self.velocities = []
        self.forces = []

    def Read(self,trajectoryFile,filename): # Needs a trajectory file reader
        return trajectoryFile.Read(self, filename)

    def Copy(self):
        import copy
        newTrj = copy.copy(self)
        newTrj.positions = []
        newTrj.velocities = []
        newTrj.forces = []
        for pos in self.positions:
            newTrj.positions.append(pos.copy())
        for vel in self.velocities:
            newTrj.velocities.append(vel.copy())
        for force in self.forces:
            newTrj.forces.append(force.copy())
        return newTrj

    def DropFrame(self,index):
        # Drop a frame by index. Useful when frames are duplicate
        try:
            del self.timesteps_of_each_frame[index]
            del self.positions[index]
            if len(self.velocities) > 0:
                del self.velocities[index]
            if len(self.forces) > 0:
                del self.forces[index]
            self.NFrames -= 1
        except:
            error("In Trajectory.DropFrame(), index {} is invalid.".format(index))
        return True




class MolecularSystem:
    # A molecular system is a collection of molecules and, in some cases, associated information including boundary
    # conditions, force field descriptions, overall thermodynamic properties, etc.
    # It can also represent the trajectory of a MD simulation if the system contains only 1 molecule, for example in
    # analyzing the results of a QM run. If the system contains multiple molecules, a collection of MolecularSystems
    # shall be needed to represent a trajectory.
    molecules: [Molecule]
    boundary: []
    origin: []
    interMolecularBonds: [Bond]
    name: str
    trajectory: Trajectory
    def __init__(self,name=None):
        self.molecules = []
        self.boundary = None    # boundary should be a 3x3 matrix for an orthogonal system
        self.origin = None      # origin should be a 3 vector
        self.interMolecularBonds = []
        self.name = name
        self.trajectory = None
    def Read(self,molecularFile,filename):
        molecularFile.Read(self,filename)
        # perform a consistency check, and create necessary auxiliary data structures
        for m in self.molecules:
            if not m.CheckConsistency():
                return False

    def ReadTrajectory(self,TrajectoryFile, filename, timestep_in_fs = 1.0,dtype='float32'):
        # timestep is set upon the first call of ReadTrajectory. Later calls do not alter the timestep
        if self.trajectory == None:
            self.trajectory = Trajectory(self,timestep_in_fs,dtype)
        self.trajectory.Read(TrajectoryFile,filename)


    def Write(self,molecularFile):
        molecularFile.Write(self)

    def Copy(self):
        # This Copy() function will deep-copy everything except the trajectory, for which only a reference is copied
        # If really want to copy everythin, use CopyWithTrajectory()
        import copy
        newMS = copy.copy(self)
        newMS.molecules = []
        newMS.interMolecularBonds = []
        for m in self.molecules:
            newMS.molecules.append(m.Copy())
        for b in self.interMolecularBonds:
            newMS.interMolecularBonds.append(b.Copy())
        newMS.boundary = copy.deepcopy(self.boundary)
        newMS.trajectory = self.trajectory   # copies only a reference
        return newMS

    def CopyWithTrajectory(self):
        newMS = self.Copy()
        newMS.trajectory = self.trajectory.Copy()
        return newMS

    def AtomCount(self):
        atomCount = 0
        for m in self.molecules:
            atomCount += len(m.atoms)
        return atomCount

    def RenumberAtomSerials(self, startingSystemwideSerial = 1):
        # renumber atom serials and systemwideSerials. This is required in most cases where the atoms must have a
        # consecutive order so that other programs can read it.
        oldToNewSerialMap = [ {} for i in range(len(self.molecules)) ]
        oldToNewSystemwideSerialMap = {}

        counterInSystem = startingSystemwideSerial
        for i,m in enumerate(self.molecules):
            counterInMolecule = 1
            for a in m.atoms:
                newSerial = '{}'.format(counterInMolecule)
                newSystemwideSerial = '{}'.format(counterInSystem)
                oldToNewSerialMap[i][a.serial] = newSerial
                if a.systemwideSerial != None:
                    oldToNewSystemwideSerialMap[a.systemwideSerial] = newSystemwideSerial
                a.serial = newSerial
                a.systemwideSerial = newSystemwideSerial
                counterInMolecule += 1
                counterInSystem += 1

        for i,m in enumerate(self.molecules):
            for b in m.bonds:
                b.atom1 = oldToNewSerialMap[i][b.atom1]
                b.atom2 = oldToNewSerialMap[i][b.atom2]

        for b in self.interMolecularBonds:
            b.atom1 = oldToNewSystemwideSerialMap[b.atom1]
            b.atom2 = oldToNewSystemwideSerialMap[b.atom2]

    def AutoDetectBonds(self, autoBondDetector, flushCurrentBonds = False):
        autoBondDetector.Detect(self, flushCurrentBonds)

    def Translate(self,dx,dy,dz):
        # Molecule.Translate is not defined. Seems unnecessary.
        for m in self.molecules:
            for a in m.atoms:
                a.Translate(dx,dy,dz)

    def FractionalToCartesianCoordinates(self):
        # Now support only orthogonal systems.
        A = self.boundary[0][0]
        B = self.boundary[1][1]
        C = self.boundary[2][2]
        for m in self.molecules:
            for a in m.atoms:
                a.x *= A
                a.y *= B
                a.z *= C

    def UpdateCoordinatesByTrajectoryFrame(self,iFrame):
        if self.trajectory == None:
            return False
        for mol in self.molecules:
            for atom in mol.atoms:
                if atom.systemwideSerial in self.trajectory.serial_to_index_map:
                    pos = self.trajectory.serial_to_index_map[atom.systemwideSerial]
                    atom.x = self.trajectory.positions[iFrame][pos][0]
                    atom.y = self.trajectory.positions[iFrame][pos][1]
                    atom.z = self.trajectory.positions[iFrame][pos][2]


    # Read and Write operations are delegated to concrete MolecularFile classes.
    def Summary(self):
        molCount = len(self.molecules)
        bondCount = 0
        for m in self.molecules:
            bondCount += len(m.bonds)
            #m.Summary()
        output("MolecularSystem {} has {} molecules, {} atoms, {} intra-molecular bonds, and {} inter-molecular bonds".format(
            self.name,molCount,self.AtomCount(),bondCount,len(self.interMolecularBonds)))
