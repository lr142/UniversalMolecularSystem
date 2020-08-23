from  UniversalMolecularSystem import *
import Utility
import math

class LAMMPSDUMPFile (MolecularFile):

    def __init__(self):
        pass

    def __ReadAFrame(self,trajectory,contents,line_no):
        # 1st line contains the timestep of current frame
        line_no += 1
        timestep = int(contents[line_no])
        if_first_frame = (True if trajectory.NFrames == 0 else False)
        trajectory.timesteps_of_each_frame.append(timestep)

        def _skip_to(flag,line_no):
            while line_no < len(contents):
                line = contents[line_no]
                if line.startswith(flag):
                    return line_no
                else:
                    line_no += 1

        line_no = _skip_to('ITEM: NUMBER OF ATOMS',line_no) + 1

        NAtomsInThisFrame = int(contents[line_no].strip())

        if if_first_frame:  # Read in NAtoms on 1st frame. Otherwise just skip to the data section
            trajectory.NAtoms = NAtomsInThisFrame
        else:
            if NAtomsInThisFrame != trajectory.NAtoms:
                error("NAtoms in frame {} with timestep {} is different from the rest!".format(trajectory.Nframes,timestep),False)
                error("{} atoms in this framehere, while others have {} atoms".format(NAtomsInThisFrame,trajectory.NAtoms),False)

        line_no = _skip_to('ITEM: ATOMS',line_no)
        line = contents[line_no]
        parts = line[12:].strip().split()
        # column number of the following properties, starting from 1
        def _index_of_keyword(key):
            return parts.index(key) if key in parts else -1

        idcol = _index_of_keyword('id')
        xcol = _index_of_keyword('x')
        ycol = _index_of_keyword('y')
        zcol = _index_of_keyword('z')
        vxcol = _index_of_keyword('vx')
        vycol = _index_of_keyword('vy')
        vzcol = _index_of_keyword('vz')
        fxcol = _index_of_keyword('fx')
        fycol = _index_of_keyword('fy')
        fzcol = _index_of_keyword('fz')

        ifVelocity = False if (vxcol == -1 or vycol == -1 or vzcol == -1) else True
        ifForce = False if (fxcol == -1 or fycol == -1 or fzcol == -1) else True

        # The dump file must at least include id and x,y,z
        if idcol == -1 or xcol == -1 or ycol == -1 or zcol == -1:
            error("The DUMP file must at least include atom id and x,y,z",False)
            return False

        # Atoms in LAMMPSDUMP file comes in disorder. Before reading the first frame, we don't know which atoms are included
        # in the dump. Therefore special care should be taken upon reading the 1st frame.


        list_ids_to_be_sorted = [-1 for _ in range(NAtomsInThisFrame)] if if_first_frame else None
        list_ids = ["" for _ in range(NAtomsInThisFrame)]
        list_xyzs = np.zeros((NAtomsInThisFrame,3))
        list_velocities = np.zeros((NAtomsInThisFrame,3))
        list_forces = np.zeros((NAtomsInThisFrame,3))



        for i in range(NAtomsInThisFrame):
            line_no += 1
            parts = contents[line_no].strip().split()
            id = parts[idcol]
            if if_first_frame:
                list_ids_to_be_sorted[i] = int(id)
            list_ids[i] = id

            list_xyzs[i][0],list_xyzs[i][1],list_xyzs[i][2] = \
                [ float(parts[xcol]), float(parts[ycol]), float(parts[zcol]) ]

            if ifVelocity:
                list_velocities[i][0], list_velocities[i][1], list_velocities[i][2] = \
                    [ float(parts[vxcol]), float(parts[vycol]), float(parts[vzcol]) ]

            if ifForce:
                list_forces[i][0], list_forces[i][1], list_forces[i][2] = \
                    [ float(parts[fxcol]), float(parts[fycol]), float(parts[fzcol]) ]

        if if_first_frame:
            list_ids_to_be_sorted.sort()
            # Reminder: the following two flags are members of Trajectory
            #     serial_to_index_map: map   # A map that maps atom.systemwideSerial to indexes in the following array
            #     index_to_serial: []        # with length NAtom, this array records each atom's systemwideSerial to its MolecularSystem
            trajectory.index_to_serial = [""] * NAtomsInThisFrame
            for index,serial in enumerate(list_ids_to_be_sorted):
                trajectory.serial_to_index_map[str(serial)] = index
                trajectory.index_to_serial[index] = str(serial)
            #Debugging
            # for i,serial in enumerate(trajectory.index_to_serial):
            #     index = trajectory.serial_to_index_map[serial]
            #     print('{}   {}    {}'.format(i,serial,index))
            #     if i != index:
            #         error("ERRORRRORRORRRORRORR")

        # Now we create the Data Structure
        array_xyzs = np.zeros((NAtomsInThisFrame,3))
        array_velocities = np.zeros((NAtomsInThisFrame,3)) if ifVelocity else None
        array_forces = np.zeros((NAtomsInThisFrame,3)) if ifForce else None


        for indexInFile, serial in enumerate(list_ids):
            position = trajectory.serial_to_index_map[serial]
            array_xyzs[position,:] = list_xyzs[indexInFile]
            if ifVelocity:
                array_velocities[position,:] = list_velocities[indexInFile]
            if ifForce:
                array_forces[position,:] = list_forces[indexInFile]

        trajectory.positions.append(array_xyzs)
        trajectory.velocities.append(array_velocities)
        trajectory.forces.append(array_forces)

        trajectory.NFrames += 1


        return line_no + 1

    def Read(self,trajectory,filename):

        if trajectory.parentMolecularSystem == None or len(trajectory.parentMolecularSystem.molecules) == 0:
            error("We don't allow a LAMMPSDUMPFile to be read alone. The caller should construct a ",False)
            error("MolecularStructure at first by reading a LAMMPS DATA file first, then call ",False)
            error("molSys.ReadTrajectory() from that MolecularStructure object instead.",False)
            return False

        file = None

        contents = None
        try:
            with open(filename,'r') as file:
                contents = file.readlines()
        except:
            error("Can't open LAMMPS DUMP file {}".format(filename),False)
            return False

        line_no = 0

        output("Reading DUMP file {}... ".format(filename))
        while True:
            if line_no >= len(contents):
                break
            elif contents[line_no].startswith('ITEM: TIMESTEP'):
                line_no = self.__ReadAFrame(trajectory,contents,line_no)
                Utility.ProgressBar(line_no/len(contents),50)
            else:
                line_no += 1
        Utility.ProgressBar(1.0, 50)
        output('')
        return True


def TestLAMMPSDUMPFile():

    from pathlib import Path
    path = '{}/{}'.format(Path.home(),'SCPDest/03.MDRuns/008.FillCH4/1.0/')

    ms = MolecularSystem()
    ms.Read(LAMMPSDATAFile(),path+'equil.data')
    ms.Summary()


    index = 0
    mask = [False] * ms.AtomCount()
    for m in ms.molecules:
        for a in m.atoms:
            if len(m.atoms) <= 5:
                mask[index] = True
            elif a.y > 75:
                mask[index] = True
            else:
                mask[index] = False
            index += 1
    subSys = ReduceSystemToOneMolecule(SubSystemByMask(ms,mask))

    ms.ReadTrajectory(LAMMPSDUMPFile(),path+'dump1.lammpstrj',timestep_in_fs=0.5)
    ms.trajectory.DropFrame(-1)
    ms.ReadTrajectory(LAMMPSDUMPFile(),path+'dump2.lammpstrj')


    output.setoutput(open('dump.xyz', 'w'))

    for iFrame in range(ms.trajectory.NFrames):
        print(iFrame)
        subSys.UpdateCoordinatesByTrajectoryFrame(iFrame)
        toWrite = ReduceSystemToOneMolecule(subSys)
        toWrite.Write(XYZFile())
