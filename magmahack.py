from subprocess import check_output,Popen, PIPE
import time
import pyperclip

from subprocess import Popen, PIPE




def Keypress(Seq):
    if type(Seq)==str:
        Seq = [Seq]
    Popen(['xdotool',"key"]+Seq)
    time.sleep(0.1)
   
def Type(Seq):
    if type(Seq)==str:
        Seq = [Seq]
    Popen(['xdotool',"type"]+Seq)
    time.sleep(.1)

def Command(Seq):
    Type(Seq)
    Keypress("Return")

def FocusWindow(id_):
    Popen(["xdotool", "windowfocus", id_])
    time.sleep(.1)

'''
SETUP:

The ~-directory on server must be mirror of home directory.
There is a terminal called "MagmaTerminal", which has a screen
session and two tabs: 
on 0, there's a magma console
on 1, there's a command line in ~
'''


def updateFile(filename="Hecke.m"):
    
    FocusWindow(WINDOWMAGMA)
    
    Keypress(["Control_L+a","1"])
    Type("vim "+filename)
    Keypress("Return")

    file = open(filename)
    text= "".join(i for i in file)
    temp = pyperclip.paste()
    pyperclip.copy(text)
    
    Type(":1, $d")
    Keypress("Return")
    
    Keypress("i")
    Keypress("Control_L+Shift+V")
    
    Keypress("Escape")
    Type(":w")
    Keypress("Return")
    Type(":q")
    Keypress("Return")
    
    Keypress(["Control_L+a","0"])
    
    pyperclip.copy(temp)
    
    FocusWindow(WINDOWNOW)
    
def setupMagma(setupscript):
    file = open(setupscript)
    text= "".join(i for i in file)
    temp = pyperclip.paste()
    pyperclip.copy(text)
    
    FocusWindow(WINDOWMAGMA)
    
    Keypress(["Control_L+a","0"])
    
    Keypress(["Control_L+x"]) #delete what's currently on command line
    Keypress("Control_L+Shift+V")
    pyperclip.copy(temp)

def setupSession():
    WINDOWMAGMA = check_output(["xdotool",'search',"--name","MagmaTerminal"])
    FocusWindow(WINDOWMAGMA)
    
    Command("screen")
    Keypress("Return")
    
    Command("cd ~/Hecke.m")
    Command("magma")
    Keypress(["Control_L+a","Control_L+c"])
    Keypress(["Control_L+a","0"])
    return

WINDOWMAGMA = check_output(["xdotool",'search',"--name","MagmaTerminal"])
WINDOWNOW = check_output(["xdotool","getwindowfocus"])
FocusWindow(WINDOWMAGMA)

updateFile("Hecke.m")
setupMagma("HeckeUpstartscript.m")
FocusWindow(WINDOWNOW)

