from MOL2File import *
from XYZFile import *
from XSDFile import *
from PDBFile import *
from BondDetection import *
from LAMMPSDATAFile import *
from LAMMPSDUMPFile import *

# sys.argv = ["","testcase/2CH4.xyz"]
# MainAsXYZ2Mol2(len(sys.argv),sys.argv)
#
#MainAsGaussianLog2XYZ(2,["","testcase/in.log"])

#TestMOL2File("testcase/Given1.mol2")

#TestXYZFile("testcase/2CH4.xyz")

#TestPDBFile()

#TestXSDFile()

#TestBondDetection()


def TestLAMMPSDUMPFile():

    ms = MolecularSystem()
    ms.Read(LAMMPSDATAFile(),'equil.data')

    ms.ReadTrajectory(LAMMPSDUMPFile(),'dump1.lammpstrj',timestep_in_fs=0.5)
    ms.ReadTrajectory(LAMMPSDUMPFile(),'dump2.lammpstrj')

    newMS = ms.CopyWithTrajectory()

    ms.trajectory.positions = []   # intentionally destroys the original structure
    ms.trajectory.NFrames = -1009

    index = 0
    print(newMS.trajectory.index_to_serial[index])
    for i in range(newMS.trajectory.NFrames):
        import sys

        sys.stdout.write('{} '.format(newMS.trajectory.positions[i][index]))

        sys.stdout.write('{} '.format(newMS.trajectory.velocities[i][index]))

        sys.stdout.write('{} '.format(newMS.trajectory.forces[i][index]))

        sys.stdout.write('\n')


    # ms.Summary()
    # output.setoutput(open('dump.mol2','w'))
    # ms.Write(MOL2File())

TestLAMMPSDUMPFile()




