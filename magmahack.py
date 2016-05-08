from subprocess import check_output,Popen, PIPE
import time
import pyperclip

from subprocess import Popen, PIPE
import sys

def Press(Seq):
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
    Press("Return")

def FocusWindow(id_):
    Popen(["xdotool", "windowfocus", id_])
    time.sleep(.1)
def RaiseWindow(id_):
    Popen(["xdotool", "windowraise", id_])

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
    
    Press(["Control_L+a","1"])
    Type("vim "+filename)
    Press("Return")

    file = open(filename)
    text= "".join(i for i in file)
    temp = pyperclip.paste()
    pyperclip.copy(text)
    
    Type(":1, $d")
    Press("Return")
    
    Press("i")
    Press("Control_L+Shift+V")
    
    Press("Escape")
    Command(":w")
    Command(":q")
    
    pyperclip.copy(temp)
    
    Press(["Control_L+a","0"])
    
    FocusWindow(WINDOWNOW)
    
def setupMagma(setupscript):
    file = open(setupscript)
    text= "".join(i for i in file)
    temp = pyperclip.paste()
    pyperclip.copy(text)
    
    FocusWindow(WINDOWMAGMA)
    
    Press(["Control_L+a","0"])
    Command(["magma -b"])  #the combination of these two is just
    time.sleep(1)
    Press(["Control_L+c"]) # to ensure there is a new magma windows
    Press(["Control_L+c"]) # with all prior calculations cancelled
    Press(["Control_L+d"]) # and then started anew
    Command(["magma -b"])
    
    Press("Control_L+Shift+V")
    pyperclip.copy(temp)

WINDOWMAGMA = check_output(["xdotool",'search',"--name","MagmaTerminal"])
WINDOWNOW = check_output(["xdotool","getwindowfocus"])
FocusWindow(WINDOWMAGMA)
RaiseWindow(WINDOWMAGMA)

file_called = sys.argv[1]
if file_called == "Hecke.m":
    updateFile(file_called)
    setupMagma("HeckeUpstartscript.m")
    FocusWindow(WINDOWNOW)
elif file_called.endswith("script"):
    updateFile(file_called)

elif file_called.startswith("Job"):
    updateFile(file_called)
    setupMagma(file_called)
else:
    setupMagma(file_called)
    
