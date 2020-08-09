from Utility import *
from UniversalMolecularSystem import *
from MolecularManipulation import *
from MOL2File import *

# A Materials Studio XSD files has a format similar to, if not simpler than, a HTML file.
# Two types of lines may be encountered:
# Type 1:
# <SectionName Property1="Value1" Property2="Value2" ... PropertyN="ValueN">
# ...... Several lines in between which provides additional information about this section
# </SectionName>
# The above line marks the end of this section
# Type 2:
# <SectionName Property1="Value1" Property2="Value2" ... PropertyN="ValueN"/>
# As you can see, in the second type the whole line is a single section. In addition, seems like Materials Studio never
# breaks a line in the middle for a pair of < and > brackets, for either type 1 and 2. This feature makes the reading of
# a XSD file easier.

# A Node is a 'Section' as mentioned above, which has a name, a bunch of properties, and 0 to many children
class XSDNode:
    name:str
    properties:map
    children: []
    def __init__(self,name):
        self.name = name
        self.properties = {}
        self.children = None

# This function creates a XSDNode from a line like below with both leading and trailing <, >, or /> symbols striped away
# SectionName Property1="Value1" Property2="Value2" ... PropertyN="ValueN". Its children information is not dealt with
# in this function. It return None if the given line has a bad format.
def NodeFromLine(line):
    try:
        result = StringTok(line,' ')

        if result == None:   # Just has a name with no properties, like <SectionName> or <SectionName/>
            name = line
            return XSDNode(name)

        (name,line) = StringTok(line,' ')

        newNode = XSDNode(name)
        while len(line) > 0:
            (key,line) = StringTok(line,'="')
            (value,line) = StringTok(line,'"')
            newNode.properties[key] = value
        return newNode
    except:
        return None

# This function reads out a molecular structure from a tree of XSDNodes. rootNode is given as the parameter.
# Returns True if successful, otherwise False
def ReadMolecularStructureFromAXSDNodeTree(molecularSystem, rootNode):
    # Will use a stack 'nodeStack' to implement a depth-first search
    nodeStack = [rootNode]
    molCount = 0
    # import note: Atoms and Molecule may not occur in a correct order. For example, a Molecule node may appear after
    # all of its constituent atoms are listed. Besides, for protein structures, atoms may belong to a Chain rather than
    # a Molecule. In the 1st pass, we first browse the nodeTree and collect all Atoms, Chains, Bonds, and Molecule. In the
    # 2nd pass, we assign each Atom to a Molecule/Chain based on the 'Parent' property of each Atom and the 'ID' of
    # of each Molecule or Chain object. In the third pass, we remove those Molecules that contain no Atoms (happens
    # when the Molecule contains only several Chains but no Atom directly).

    tempMoleculesCollection = []
    tempAtomsCollection = []
    tempBondCollection = []
    # 1st pass
    while len(nodeStack) > 0:
        node = nodeStack.pop()
        if node.name == 'Molecule' or node.name == 'Chain' or node.name == 'SubUnit':
            newMolecule = Molecule()
            newMolecule.name = node.properties['Name']
            newMolecule.serial = node.properties['ID']
            tempMoleculesCollection.append(newMolecule)
        elif node.name == 'Atom3d':
            newAtom = Atom()
            newAtom.element = node.properties['Components']
            newAtom.serial = newAtom.systemwideSerial = node.properties['ID']
            (newAtom.x,newAtom.y,newAtom.z) = (float(num) for num in node.properties['XYZ'].split(','))

            if 'Parent' in node.properties:
                newAtom.parent = node.properties['Parent']
            else:  # If the Atom doesn't have a 'Parent' property, the lastly added molecule is its parent
                newAtom.parent = tempMoleculesCollection[-1].serial


            # Name and Hybridization properties are not always present
            newAtom.name = node.properties['Name'] if 'Name' in node.properties else '{}{}'.format(newAtom.element,newAtom.serial)
            newAtom.type = "{}{}".format(newAtom.element,'.'+node.properties['Hybridization'] if 'Hybridization' in node.properties else '')
            tempAtomsCollection.append(newAtom)
        elif node.name == 'Bond':
            newBond = Bond()
            (newBond.atom1,newBond.atom2) = node.properties['Connects'].split(',')
            newBond.type = node.properties['Type'] if 'Type' in node.properties else '1'    # Type is also not always present
            tempBondCollection.append(newBond)

        elif node.name == 'SpaceGroup':   # Unit Cell size can be read. The info are given like below. We need the AVector/BVector/CVector
            #<SpaceGroup ID="4" Parent="2" Children="83672" AVector="79.1697982956743,0,0" BVector="0,79.1697982956743,0" CVector="0,0,79.1697982956743"
            # OrientationBase="C along Z, B in YZ plane" Centering="3D Primitive-Centered"
            # Lattice="3D Triclinic" GroupName="P1" Operators="1,0,0,0,0,1,0,0,0,0,1,0" DisplayRange="0,1,0,1,0,1" LineThickness="2"
            # CylinderRadius="0.2" LabelAxes="1" ActiveSystem="2" ITNumber="1" LongName="P 1" Qualifier="Origin-1" SchoenfliesName="C1-1"
            # System="Triclinic" Class="1"/>
            molecularSystem.boundary = [[],[],[]]
            vectors = ['AVector','BVector','CVector']
            for i, vec in enumerate(vectors):
                v = node.properties[vec].split(',')
                for j in range(3):
                    v[j] = float(v[j])
                import copy
                molecularSystem.boundary[i] = copy.copy(v)

        else:
            pass

        if node.children != None:
            node.children.reverse()
            nodeStack.extend(node.children)  # Must add children reversely to stack in order to preserve the correct order

    # 2nd pass, assign atoms and bonds to each molecule
    moleculeSerialToIndexMap = {}
    for i,m in enumerate(tempMoleculesCollection):
        moleculeSerialToIndexMap[m.serial] = i

    serialToAtomMap = {}
    for atom in tempAtomsCollection:
        serialToAtomMap[atom.serial] = atom
        molID = atom.parent;
        molecule = tempMoleculesCollection[moleculeSerialToIndexMap[molID]]
        molecule.atoms.append(atom)

    for bond in tempBondCollection:
        atom1 = serialToAtomMap[bond.atom1]
        atom2 = serialToAtomMap[bond.atom2]
        if atom1.parent == atom2.parent:      # intra-molecular bond
            molID = atom1.parent
            tempMoleculesCollection[moleculeSerialToIndexMap[molID]].bonds.append(bond)
        else:   # inter-molecular bond
            molecularSystem.interMolecularBonds.append(bond)

    # 3rd pass, remove (forget) those empty molecules
    for mol in tempMoleculesCollection:
        if len(mol.atoms) > 0:
            molecularSystem.molecules.append(mol)

    # In addition, we must renumber the atom serials, so that other programs can read the structure
    molecularSystem.RenumberAtomSerials()


    return True

class XSDFile(MolecularFile):

    def __init__(self):
        pass

    def Read(self,molecularSystem,filename):
        contents = None
        try:
            with open(filename,'r') as file:
                contents = file.readlines()
        except:
            error("Can't open XSD file {}".format(filename),False)
            return False

        rootNode = XSDNode("root")
        rootNode.children = []
        nodeStack = []
        nodeStack.append(rootNode)  # This list is used as a stack, which records the top level element

        for i,line in enumerate(contents):
            line = line.strip()
            if len(line) == 0:
                continue

            newNode = None

            if line.startswith("<?") or line.startswith("<!"):   # Special line at the beginning. Simply ignore them.
                continue

            elif line.startswith('</'):  # End section of Type 1
                nodeStack.pop()  # pop from the Stack
                continue

            elif line.startswith('<') and line.endswith('/>'):  # Type 2
                newNode = NodeFromLine(line[1:-2])
                if newNode != None:
                    newNode.children = None
                    nodeStack[-1].children.append(newNode)

            elif line.startswith('<') and line.endswith('>'):  # Type 1
                newNode = NodeFromLine(line[1:-1])
                # Tip: When a notation like string[A:-B] is used, it simply means 'cut off A and B characters from
                # the beginning and the end'
                if newNode != None:
                    newNode.children = []
                    # Push into the stack.
                    nodeStack[-1].children.append(newNode)
                    nodeStack.append(newNode)
            else:
                pass

            # Check error at this round:
            if newNode == None:
                error('Unexpected format @ line {} of XSD file {}:\n{}'.format(i+1,filename,line))

        # DFSSearch is for debugging purposes: show the tree structure of what just read
        def DFSSearch(node, level,hierarchy):
            blanks = ''.join( ['  ' for i in range(level)])
            output("{}HIERARCHY: {}".format(blanks,hierarchy))
            output("{}NAME: {}".format(blanks,node.name))
            propertyString = "{}PROPERTIES: {}".format(blanks, "None" if len(node.properties) == 0 else "")
            for p in node.properties:
                propertyString += "{}=\"{}\" ".format(p,node.properties[p])
            output(propertyString)
            if node.children!=None:
                output("{}CHILDREN:".format(blanks))
                for i,c in enumerate(node.children,start=1):
                    DFSSearch(c,level+1,"{}-{}".format(hierarchy,i))

        # Seems like the hierarchy is always ROOT -> XSD -> AtomisticTreeRoot
        # Therefore to simply it a little bit, we start searching from this node (AtomisticTreeRoot)
        # output.setoutput(open('dump.txt','w'))
        # DFSSearch(rootNode.children[0].children[0],0,1)

        return ReadMolecularStructureFromAXSDNodeTree(molecularSystem,rootNode)



    def Write(self,molecularSystem):
        error("Sorry! Writing in the XSD format is not supported!")

def TestXSDFile():
    ms = MolecularSystem()
    ms.Read(XSDFile(),'testcase/Coal100.xsd')

    ms.Summary()
    output(ms.boundary)

    output.setoutput(open('testcase/dump.mol2', 'w'))
    ms.Write(MOL2File())

