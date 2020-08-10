# Public utility functions and classes. Implement commonly used constants, mathematical functions, etc
import math
import sys

class PeriodicTable:
    def __init__(self):
        self.elements = "M H He "\
        "Li Be B C N O F Ne "\
        "Na Mg Al Si P S Cl Ar "\
        "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "\
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "\
        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "\
        "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus Uuo";
        self.elements = self.elements.split()
    # M is the dummy atom

    def AtomicNumberToElement(self,i):
        if i >= len(self.elements):
            i = 0
        return self.elements[i]

    def ElementToAtomicNumber(self,element):
        for i in range(len(self.elements)):
            if element == self.elements[i]:
                return i
        return -1

class ErrorHandler():
    # The ErrorHandling class controls how to output and deal with an error message
    # There exists only ONE ErrorHandling object which is called 'error'. Here we used the 'singleton' design pattern
    # The 'error' object must be constructed once and only once
    def __init__(self):
        self.outputfile = None
        self.OnOffState = True

    def setoutput(self,outputfile):
        self.outputfile = outputfile

    def __call__(self,message,fatal = True):

        # Warnings can be turned off. Fatal errors can't be turned off
        if (not self.OnOffState) and (not fatal):
            return

        message = "Fatal Error: "+message if fatal else "Warning: "+message
        print(message) # output to stdout
        if self.outputfile != None: # and to the outputfile. If no output file is set up, to stdout only
            print(message, file = self.outputfile)
        if fatal:
            quit()

    def turn_on(self):
        self.OnOffState = True

    def turn_off(self):
        self.OnOffState = False

error = ErrorHandler()  # A global error-handling object (the only object of this class).
# in main.py, error.setoutput() can be called to direct the error message into either the command line screen or somewhere
# else (for example in a GUI version, direct the error message to a message box)

class OutputHandler():
    def __init__(self):
        self.outputfile = None
        self.OnOffState = True

    def setoutput(self, outputfile):
        if self.outputfile != None:
            self.outputfile.close()
        self.outputfile = outputfile

    def __call__(self,message):

        if self.OnOffState == False: # Output is turned off
            return

        if self.outputfile != None:
            print(message,file = self.outputfile)
        else:
            print(message)

    def turn_on(self):
        self.OnOffState = True

    def turn_off(self):
        self.OnOffState = False

output = OutputHandler() # A unique and global object, similar as the ErrorHandler()

# Mathematical Functions
def Distance(atom1, atom2):
    return math.sqrt(
        math.pow(atom1.x - atom2.x,2) +
        math.pow(atom1.y - atom2.y,2) +
        math.pow(atom1.z - atom2.z,2)
    )


# String Manipulations
# Note: Due to my lack of knowledge in this field, there may be better ways to implement the following functions,
# for example through regular expressions. Please bear with the clumsiness of these functions. At least they work.

def StringTok(string, token):
    # Finds the first occurrence of the token in the string and split the string into a tuple (first,second) without
    # the given token. For example StringTok("xyz = 543",'=') returns ("xyz" and "543"). (Both parts are with
    # beginning and ending blank spaces striped away. If no such token is found, it returns None.
    pos = string.find(token)
    if pos == -1:
        return None
    else:
        return (string[0:pos].strip(), string[pos + len(token):].strip())


def ProgressBar(percent,length = 50):
    filled = round(length * percent)
    bar = '['+''.join(['#' for i in range(filled)]) + ''.join([' ' for i in range(length-filled)])+']'
    sys.stdout.write('\r{} {:.2f}%'.format(bar,percent*100.0))
