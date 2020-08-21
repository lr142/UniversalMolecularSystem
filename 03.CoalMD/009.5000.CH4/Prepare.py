import sys
sys.path.append("../")
import StructurePrepare

para = StructurePrepare.StructruePrepParameter()

para.numberOfCoalCellsInXDirection = 2
para.extraPaddingInXDirection = 0
para.outerDiameter = 150.0
para.innerDiameter = 100.0
para.lengthOfFilledCH4 = 198  # fill the coal tube with CH4 in range [0, lengthOfFilledCH4)
para.lengthOfFilledH2O = 0


para.methaneCellFile = '{}/{}/{}'.format("../000.Common.Files.Read.Only", "water.and.ch4.cells", "5000methane.pdb")
# Set a loose criteria for methane
para.minDistanceForClashing = 3.0

StructurePrepare.Prepare(para)