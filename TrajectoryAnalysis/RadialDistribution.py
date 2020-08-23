
def Read():

    from pathlib import Path
    path = '{}/{}'.format(Path.home(),'SCPDest/03.MDRuns/008.FillCH4/1.0/')

    ms = MolecularSystem()
    ms.Read(LAMMPSDATAFile(),path+'equil.data')

    ms.ReadTrajectory(LAMMPSDUMPFile(),path+'dump1.lammpstrj',timestep_in_fs=0.5)
    ms.trajectory.DropFrame(-1)
    ms.ReadTrajectory(LAMMPSDUMPFile(),path+'dump2.lammpstrj')

    return ms


