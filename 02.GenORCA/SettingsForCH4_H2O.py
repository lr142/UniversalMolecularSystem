xyzfile = "ch4+h2o.xyz"
# CommandLine = "B3LYP D3 def2-TZVP def2/J RIJCOSX noautostart miniprint nopop"
# CommandLine = "PWPB95 D3 def2-QZVPP def2/J def2-QZVPP/C RIJCOSX grid4 gridx4 tightSCF noautostart miniprint nopop"
CommandLine = "DLPNO-CCSD(T) normalPNO RIJK cc-pVTZ cc-pVTZ/JK cc-pVTZ/C tightSCF noautostart miniprint nopop"
maxcore = 3500
nprocs = 10
part1_atoms = set(range(5))  # atom indexes starting from 0
part2_atoms = set(range(5, 8))

distance = list(range(20, 81, 2))  # in units of 0.1 Å
for i, d in enumerate(distance):
    distance[i] /= 10.0

distance_measure_part1 = [0]
# atom indexes in each part for measuring distances. On part1, multiple atoms may be used as the indicator, which
# is useful in measuring the distance between a water and an aromatic ring.
distance_measure_part2 = 5

prefix = "ch4h2o"