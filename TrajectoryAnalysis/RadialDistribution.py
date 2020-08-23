import sys
sys.path.append('../')
from UniversalMolecularSystem import *
from LAMMPSDUMPFile import LAMMPSDUMPFile
from MOL2File import MOL2File
from XYZFile import XYZFile
from Utility import *
import MolecularManipulation
import BondDetection


def Read():

    from pathlib import Path
    path = '{}/{}'.format(Path.home(),'SCPDest/03.MDRuns/008.FillCH4/1.0/')

    ms = MolecularSystem()
    ms.Read(MOL2File(),path+'../positions.mol2')
    ms.Summary()

    ms.ReadTrajectory(LAMMPSDUMPFile(),path+'dump1.lammpstrj',timestep_in_fs=0.5)
    # ms.trajectory.DropFrame(-1)
    # ms.ReadTrajectory(LAMMPSDUMPFile(),path+'dump2.lammpstrj')

    return ms

def SelectOnlyMethane(ms):
    mask = [False] * ms.AtomCount()
    index = 0
    for m in ms.molecules:
        for a in m.atoms:
            if len(m.atoms) > 5:
                mask[index] = False
            else:
                mask[index] = True
            index += 1
    return MolecularManipulation.SubSystemByMask(ms,mask)

def SelectOnlyMethaneCarbon(ms):
    mask = [False] * ms.AtomCount()
    index = 0
    for m in ms.molecules:
        for a in m.atoms:
            if len(m.atoms) == 5 and a.element == 'C':
                mask[index] = True
            index += 1
    return MolecularManipulation.SubSystemByMask(ms,mask)

def SelectFurtherHalfCoalTube(ms):
    mask = [False] * ms.AtomCount()
    index = 0
    for m in ms.molecules:
        for a in m.atoms:
            if len(m.atoms) <= 5:
                mask[index] = False
            elif a.y >= 75:
                mask[index] = True
            else:
                mask[index] = False
            index += 1
    return MolecularManipulation.SubSystemByMask(ms,mask)

def RadialDistribution(ms_with_trj, length, centerY, centerZ, radius, binSize, frames):
    Nbins = int(math.floor(radius/binSize))
    bins = [0.0 for _ in range(Nbins)]
    volumes = [0 for _ in range(Nbins)]   # the volume of shells
    for ibin in range(Nbins):
        rout = radius - ibin * binSize
        rin = rout - binSize if ibin < Nbins else 0.0
        volumes[ibin] = math.pi * length * ( rout * rout - rin * rin ) * 1000.0 # in unit of nm^3

    for iFrame in frames:
        ms_with_trj.UpdateCoordinatesByTrajectoryFrame(iFrame)
        for mol in ms_with_trj.molecules:
            for a in mol.atoms:
                dist_to_axis = math.sqrt(math.pow(a.y-centerY,2) + math.pow(a.z-centerZ,2))
                dist_to_edge = radius - dist_to_axis
                which_bin = int(math.floor(dist_to_edge/binSize))
                if which_bin >= binSize:
                    which_bin = -1.0
                bins[which_bin] += 1.0

    # number in each bin are averaged by NFrames
    for ibin in range(bins):
        bins[ibin] /= len(frames)

    # And further converted to NAtoms/nm^3
    for ibin in range(bins):
        bins[ibin] = bins[ibin]/volumes[ibin]
        

        kalksjdf;aj;klsjdfklajksdjfkalsdf Debug Here








def Main():
    ms_whole = Read()

    def _write_trajectory():
        ms_half_tube = SelectFurtherHalfCoalTube(ms_whole)
        ms_half_tube.Summary()
        ms_methane = SelectOnlyMethane(ms_whole)
        ms_methane.Summary()
        file = open('dump.xyz','w')
        output.setoutput(file)
        for iFrame in range(0,ms_whole.trajectory.NFrames,5):
            ProgressBar(iFrame*1.0/ms_whole.trajectory.NFrames)
            ms_methane.UpdateCoordinatesByTrajectoryFrame(iFrame)
            to_write = ms_half_tube.Copy()
            MolecularManipulation.ExtendSystem(to_write,ms_methane)
            MolecularManipulation.ReduceSystemToOneMolecule(ms_methane).Write(XYZFile())
        ProgressBar(1.00);output('')
        file.close()
    #_write_trajectory()

    ms_C = MolecularManipulation.SubSystemByMask(ms_whole)

    RadialDistribution(ms_C,centerY=75.0,centerZ=75.0,radius=53.0,binSize=0.2,frames=[-1])




if __name__ == '__main__':
    Main()


