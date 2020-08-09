from Utility import *
from XYZFile import *
from UniversalMolecularSystem import *
from MolecularManipulation import *

class PDBFile(MolecularFile):
    def __init__(self):
        pass

    # In a PDB file there may be more than one MODEL that describes the same system with different coordinates.
    # Here we just implement a simplified version that considers only one MODEL
    # We partition the system into 'Molecule's according the residue seqs. The chain information is stored in the 'parent'
    # property of each atom.
    # Bonds may be given in the 'CONECT' (not a misspell) section, or may be absent. In most cases, only bonds that are
    # not in the peptide are given. For simplicity, we just ignore all bonds and use auto-detection to find bonds.

    def Read(self,molecularSystem,filename):
        contents = None
        try:
            with open(filename,"r") as file:
                contents = file.readlines()
        except:
            error("Can't open PDB file <{}> to read.".format(filename),False)
            return False

        tempAtomCollection = []
        firstAppearedAltLoc = None
        # 1st pass, read atom info from ATOM and HETATM sections
        try:
            for i,line in enumerate(contents,start=1):
                # Pad the line with a leading blank so that reading is easier:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # The following format is found for the ATOM/HETATM section specification (version 3.30)
                    # 1 -  6      Record name   "ATOM  "
                    # 7 - 11      Integer       serial       Atom  serial number.
                    # 13 - 16      Atom      name       Atom name.
                    # 17         Character     altLoc       Alternate location indicator.
                    # 18 - 20      Residue name  resName      Residue name.
                    # 22         Character     chainID      Chain identifier.
                    # 23 - 26      Integer       resSeq       Residue sequence number.
                    # 27         AChar       iCode      Code for insertion of residues.
                    # 31 - 38      Real(8.3)     x        Orthogonal coordinates for X in Angstroms.
                    # 39 - 46      Real(8.3)     y        Orthogonal coordinates for Y in Angstroms.
                    # 47 - 54      Real(8.3)     z        Orthogonal coordinates for Z in Angstroms.
                    # 55 - 60      Real(6.2)     occupancy    Occupancy.
                    # 61 - 66      Real(6.2)     tempFactor   Temperature  factor.
                    # 77 - 78      LString(2)    element      Element symbol, right-justified.
                    # 79 - 80      LString(2)    charge       Charge  on the atom.
                    line = ' {}'.format(line.strip())
                    altLoc = line[17]  # If there is altLoc flag, only the first altLocation, such as 'A', shall be read
                    if altLoc != ' ':
                        if firstAppearedAltLoc == None:
                            firstAppearedAltLoc = altLoc
                        if altLoc != firstAppearedAltLoc:
                            continue

                    atom = Atom()
                    atom.serial = line[7:12].strip()
                    atom.name = line[13:17].strip()
                    resName = line[18:21].strip()
                    chainID = line[22]
                    resSeq = line[23:27].strip()
                    # temporarily store the resName, chainID, and resSeq in atom.parent. They will be processed later.
                    atom.parent = '{},{},{}'.format(resName,chainID,resSeq)
                    atom.x = float(line[31:39])
                    atom.y = float(line[39:47])
                    atom.z = float(line[47:55])
                    atom.element = line[77:79].strip()
                    atom.type = atom.element
                    try:
                        atom.charge = float(line[79:81])   # Charge info may be absent
                    except:
                        atom.charge = 0.0
                    tempAtomCollection.append(atom)

                else:
                    continue
        except:
            error("Unexpected format @ line {} of file {}:\n{}".format(i,filename,line),False)
            return False

        # 2nd pass, create molecules based on residue sequences and serials
        tempMoleculesCollection = {}    # Use a map to store the molecules.
        # 'ChainName+ResidueSeq' is used as the key to avoid duplication.
        for atom in tempAtomCollection:
            (resName,chainID,resSeq) = atom.parent.split(',')
            # if there is no chain info in the file, chainID here shall be ' '. This won't affect the program
            molSerial = '{}{}'.format(chainID,resSeq)
            if molSerial in tempMoleculesCollection:
                tempMoleculesCollection[molSerial].atoms.append(atom)
            else:
                newMolecule = Molecule()
                newMolecule.serial = molSerial
                newMolecule.name = resName
                tempMoleculesCollection[molSerial] = newMolecule
                newMolecule.atoms.append(atom)
            # Now we leave only the chainID in each atom. This is not an elegant treatment, maybe revise it later.
            atom.parent = chainID

        for m in tempMoleculesCollection:
            molecularSystem.molecules.append(tempMoleculesCollection[m])

        molecularSystem.RenumberAtomSerials()

        return True

    def Write(self):
        # Writing in PDB format is also not yet implemented
        error("Sorry, writing in PDB format is not supported yet!", False)
        return False


def TestPDBFile():
    ms = MolecularSystem()
    ms.Read(PDBFile(),'testcase/6zju.pdb')
    ms.Summary()
    newMS = MolecularSystem()
    newMS.molecules.append(ReduceSystemToOneMolecule(ms))
    newMS.Summary()

    ms = newMS
    output.setoutput(open('testcase/dump.mol2','w'))
    ms.Write(MOL2File())
