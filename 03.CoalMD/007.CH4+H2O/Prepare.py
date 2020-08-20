import sys
sys.path.append("../")
import StructurePrepare

para = StructurePrepare.StructruePrepParameter()

para.numberOfCoalCellsInXDirection = 2
para.extraPaddingInXDirection = 0
para.outerDiameter = 100.0
para.innerDiameter = 50.0
para.lengthOfFilledCH4 = 98  # fill the coal tube with CH4 in range [0, lengthOfFilledCH4)
para.lengthOfFilledH2O = 98

StructurePrepare.Prepare(para)