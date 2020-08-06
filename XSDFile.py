from Utility import *
from UniversalMolecularSystem import *

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

class XSDNode:
    name:str
    properties:map
    children: []
    def __init__(self,name):
        self.name = name
        self.properties = {}
        self.children = None

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

        def NodeFromLine(line):
            parts = line.split()
            name = parts[0]
            newNode = XSDNode(name)
            for index in range(1, len(parts)):
                key_and_value = parts[index].split('=')
                key = key_and_value[0]
                value = key_and_value[1][1:-1]  # [1:-1] skips the " at the beginning and the end
                newNode.properties[key] = value
            return newNode

        for i,line in enumerate(contents):
            line = line.strip()
            if len(line) == 0:
                continue

            if len(line) < 4 or line[0] != '<':
                error('Unexpected format @ line {} of XSD file {}:\n{}'.format(i + 1, filename, line),False)
                return False

            if line[0] == '<' and line[1] != '/' and line[-2] != '/': # Type 1
                newNode = NodeFromLine(line[1:-1])
                # Tip: When a notation like string[A:-B] is used, it simply means 'cut off A and B characters from
                # the beginning and the end'
                newNode.children = []
                # Push into the stack.
                nodeStack[-1].children.append(newNode)
                nodeStack.append(newNode)

            elif line[0] == '<' and line[1] != '/' and line[-2] == '/': # Type 2
                newNode = NodeFromLine(line[1:-2])
                newNode.children = None
                nodeStack[-1].children.append(newNode)
            elif line[0] == '<' and line[1] == '/': # End section of Type 1
                nodeStack.pop()  # pop from the Stack
            else:
                error('Unexpected format @ line {} of XSD file {}:\n{}'.format(i+1,filename,line))

        # Debugging: show the tree structure of what just read
        







def TestXSDFile():
    ms = MolecularSystem()
    ms.Read(XSDFile(),'testcase/Given1.xsd')