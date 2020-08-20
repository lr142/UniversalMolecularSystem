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

StructurePrepare.Prepare(para)