from MOL2File import *
from XYZFile import *

# sys.argv = ["","testcase/2CH4.xyz"]
# MainAsXYZ2Mol2(len(sys.argv),sys.argv)
#
#MainAsGaussianLog2XYZ(2,["","testcase/in.log"])

#TestMOL2File("testcase/Given1.mol2")
#TestXYZFile("testcase/2CH4.xyz")

a = Atom()
a.name = "old"
a.info = [1,2,3]

import copy
b = copy.deepcopy(a)
b.name = "new"
b.info.append(4)

print(a.name)
print(a.info)

print(b == a)
b = a
print(b == a)
