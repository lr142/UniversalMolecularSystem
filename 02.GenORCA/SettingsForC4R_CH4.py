xyzfile = "c4r+ch4.xyz"
# CommandLine = "B3LYP D3 def2-TZVP def2/J RIJCOSX noautostart miniprint nopop"
# CommandLine = "PWPB95 D3 def2-QZVPP def2/J def2-QZVPP/C RIJCOSX grid4 gridx4 tightSCF noautostart miniprint nopop"
CommandLine = "DLPNO-CCSD(T) normalPNO RIJK cc-pVTZ cc-pVTZ/JK cc-pVTZ/C tightSCF noautostart miniprint nopop"
maxcore = 3500
nprocs = 10
part1_atoms = set(range(26))  # atom indexes starting from 0
part2_atoms = set(range(26, 31))

distance = list(range(30, 61, 2))  # in units of 0.1 Å
for i, d in enumerate(distance):
    distance[i] /= 10.0

distance_measure_part1 = [2, 3, 9, 10]
# atom indexes in each part for measuring distances. On part1, multiple atoms may be used as the indicator, which
# is useful in measuring the distance between a water and an aromatic ring.
distance_measure_part2 = 26

prefix = "c4rch4"