from UniversalMolecularSystem import *
from XYZFile import *

def CalculatePlacement(ms,part2_atoms,distance,distance_measure_part1,distance_measure_part2):
    from_point = [0,0,0]
    for i in distance_measure_part1:
        from_point[0] += ms.molecules[0].atoms[i].x
        from_point[1] += ms.molecules[0].atoms[i].y
        from_point[2] += ms.molecules[0].atoms[i].z
    for i in range(3):
        from_point[i] /= len(distance_measure_part1)

    to_point = [0,0,0]
    to_atom = ms.molecules[0].atoms[distance_measure_part2]
    to_point[0] = to_atom.x
    to_point[1] = to_atom.y
    to_point[2] = to_atom.z

    import math
    def Dist(point1,point2):
        sum = 0.0
        for i in range(3):
            sum += pow(point2[i]-point1[i],2)
        return math.sqrt(sum)

    original_distance = Dist(from_point,to_point)

    output("original_distance = {}".format(original_distance))

    translation = [0,0,0]
    for i in range(3):
        translation[i] = from_point[i] + (to_point[i]-from_point[i])*distance/original_distance - to_point[i]

    for a in part2_atoms:
        ms.molecules[0].atoms[a].x += translation[0]
        ms.molecules[0].atoms[a].y += translation[1]
        ms.molecules[0].atoms[a].z += translation[2]

    newDist = Dist(from_point,[to_atom.x,to_atom.y,to_atom.z])
    output("  new_distance    = {}".format(newDist))



def Main():
    xyzfile = "fs.xyz"
    #CommandLine = "B3LYP D3 def2-TZVP def2/J RIJCOSX noautostart miniprint nopop"
    #CommandLine = "PWPB95 D3 def2-QZVPP def2/J def2-QZVPP/C RIJCOSX grid4 gridx4 tightSCF noautostart miniprint nopop"
    CommandLine = "DLPNO-CCSD(T) normalPNO RIJK cc-pVTZ cc-pVTZ/JK cc-pVTZ/C tightSCF noautostart miniprint nopop"
    maxcore = 3500
    nprocs = 12
    part1_atoms = set(range(26))    # atom indexes starting from 0
    part2_atoms = set(range(26,31))

    distance = list(range(20,61,2))   # in units of 0.1 Å
    for i,d in enumerate(distance):
        distance[i] /= 10.0

    #distance.append(3.47716)

    distance_measure_part1 = [2, 3, 9, 10]
    # atom indexes in each part for measuring distances. On part1, multiple atoms may be used as the indicator, which
    # is useful in measuring the distance between a water and an aromatic ring.
    distance_measure_part2 = 26


    ms = MolecularSystem()
    ms.Read(XYZFile(),xyzfile)

    import os
    os.system("if [ -e combined.xyz ]; then rm combined.xyz;fi")

    for dist in distance:

        newMS = ms.Copy()
        CalculatePlacement(newMS,part2_atoms,dist,distance_measure_part1,distance_measure_part2)

        lines = []
        lines.append("%pal nprocs  {} end".format(nprocs))
        lines.append("")

        for round in range(3):
            lines.append("! {}{}".format(CommandLine,"" if round==0 else " Pmodel"))
            lines.append("%maxcore {}".format(maxcore))
            lines.append("* xyz   0  1")
            for i,a in enumerate(newMS.molecules[0].atoms):
                flag = None
                if round == 1 and i in part1_atoms:
                    flag = ":"
                elif round == 2 and i in part2_atoms:
                    flag = ":"
                else:
                    flag = " "
                lines.append("{}{}       {}      {}       {}".format(a.element,flag,a.x,a.y,a.z))
            lines.append(" *")
            if round != 2:
                lines.append("\n$new_job")

            import os
            os.system('if [ ! -d {} ];then mkdir {};fi'.format(dist,dist))
            with open('{}/{}.inp'.format(dist,xyzfile.split('.')[0]),"w") as outputfile:
                for l in lines:
                    outputfile.write(l+"\n")

        currentXYZ = "{}/{}".format(dist,xyzfile)
        with open(currentXYZ,'w') as xyzhandle:
            output.setoutput(xyzhandle)
            newMS.Write(XYZFile())
        output.setoutput(None)

        os.system("cat {} >> combined.xyz".format(currentXYZ))


if __name__ == '__main__':
    Main()