from Utility import *
from UniversalMolecularSystem import *
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
    while len(nodeStack) > 0:
        node = nodeStack.pop()
        if node.name == 'Molecule':
            molCount += 1
            newMolecule = Molecule()
            newMolecule.name = node.properties['Name']
            newMolecule.serial = '{}'.format(molCount)   # You wont' get a serial for molecules
            # in an XSD file. So just number them sequentially.
            molecularSystem.molecules.append(Molecule())
        elif node.name == 'Chain':  # There is an alternative situation that atoms are not a molecule but in several
            # 'Chains', and each 'Chain' is a child of the whole 'Molecule'. In this case, we treat each chain as a
            # Molecule. This is usually true when reading a protein structure.
            newMolecule = Molecule()
            newMolecule.name = node.properties['Name']
            # Now we perform a check. If the previously added 'Molecule' has zero atoms, then all atoms should go to a
            # Chain, therefore there's no point to keep the Molecule as a level of hierarchy. We just use the Chain to
            # replace the top level Molecule. Otherwise, there may be atoms that are in the system but don't belong to
            # any chain, so we keep the top level Molecule as a separate collection of atoms.
            if len(molecularSystem.molecules[-1].atoms) == 0:
                newMolecule.serial = '{}'.format(molCount)
                molecularSystem.molecules[-1] = newMolecule;
            else:
                molCount += 1
                newMolecule.serial = '{}'.format(molCount)
                molecularSystem.molecules.append(newMolecule)
        elif node.name == 'Atom3d':
            newAtom = Atom()
            newAtom.element = node.properties['Components']
            newAtom.serial = node.properties['ID']
            (newAtom.x,newAtom.y,newAtom.z) = (float(num) for num in node.properties['XYZ'].split(','))

            # Name and Hybridization properties are not always present
            newAtom.name = node.properties['Name'] if 'Name' in node.properties else '{}{}'.format(newAtom.element,newAtom.serial)
            newAtom.type = "{}{}".format(newAtom.element,'.'+node.properties['Hybridization'] if 'Hybridization' in node.properties else '')
            molecularSystem.molecules[-1].atoms.append(newAtom)
        elif node.name == 'Bond':
            newBond = Bond()
            (newBond.atom1,newBond.atom2) = node.properties['Connects'].split(',')
            newBond.type = node.properties['Type'] if 'Type' in node.properties else '1'    # Type is also not always present
            molecularSystem.molecules[-1].bonds.append(newBond)
        else:
            pass

        if node.children != None:
            node.children.reverse()
            nodeStack.extend(node.children)  # Must add children reversely to stack in order to preserve the correct order

    # Finally, we must renumber the atom serials, so that other programs can read the structure
    oldToNewSerialMap = {}
    oldToNewSystemwideSerialMap = {}
    counterInSystem = 1
    for m in molecularSystem.molecules:
        counterInMolecule = 1
        for a in m.atoms:
            newSerial = '{}'.format(counterInMolecule)
            newSystemwideSerial = '{}'.format(counterInSystem)
            oldToNewSerialMap[a.serial] = newSerial
            oldToNewSystemwideSerialMap[a.serial] = newSystemwideSerial
            a.serial = newSerial
            a.systemwideSerial = newSystemwideSerial
            counterInMolecule += 1
            counterInSystem += 1

    for m in molecularSystem.molecules:
        for b in m.bonds:
            try:
                b.atom1 = oldToNewSerialMap[b.atom1]
                b.atom2 = oldToNewSerialMap[b.atom2]
            except:
                error("In ReadMolecularStructureFromAXSDNodeTree(): Some local bond becomes global bond")

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

    output.setoutput(open('dump.mol2','w'))
    ms.Write(MOL2File())
