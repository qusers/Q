import sys
import subprocess
import os
import glob
import shutil

def query_yes_no(question, default="no"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
            
def init():
    print("================================\n")
    print("This is the installer for Q-GPU\n")
    print("================================\n")
    print("Checking Python version")

def check_python():
    print("================================\n")
    print("This is the installer for Q-GPU\n")
    print("================================\n")
    print("Checking Python version")
    print("================================")
    print("version is {}.{}.{}".format(sys.version_info[0],
                                       sys.version_info[1],
                                       sys.version_info[2],
                                      ))
    print("================================")
    
    if sys.version_info[0] < 3: 
        print("Q-GPU only works with Python 3.6 or higher")
        print("Exiting now")
        sys.exit()

    if sys.version_info[0] == 3 and sys.version_info[1] < 6:
        print("Q-GPU only works with Python 3.6 or higher")
        print("Exiting now")    
        sys.exit()
        
def check_CUDA():
    print("================================")
    print("Checking CUDA version")
    out = subprocess.check_output(["which", "nvcc"]).decode("utf-8")

    if len(out) == 0:
        print("Can't find nvcc, exiting now")
        print("Exiting now")    
        sys.exit()

    out = subprocess.check_output(["nvcc", "--version"]).decode("utf-8")
    print("version is {}".format(out.split()[-1]))

    print("================================")
    
def get_installdir():
    print("Please provide the installation directory: ")
    INSTALLDIR = input()
    
    return INSTALLDIR

def create_installdir(INSTALLDIR):
    if os.path.isdir(INSTALLDIR) == True:
        write = query_yes_no("Directory exists, are you sure you want to install here?")
        
        if write == False:
            print("Will not write in specified folder, exiting now")
            sys.exit()
            
        write_in_folder = os.access(INSTALLDIR, os.W_OK)
    
        if write_in_folder == False or write == False:
            print("Can't write in specified folder, exiting now")
            sys.exit()
        
    else:
        try:
            os.mkdir(INSTALLDIR)
            
        except:
            print("Can't create installation directory, exiting now")
            sys.exit()
    
    ROOTDIR = os.path.abspath(INSTALLDIR)
    return ROOTDIR

def install(ROOTDIR):
    print("Installing Q-GPU in {}".format(ROOTDIR))
    for file in glob.glob('*'):
        if file == 'install.py':
            continue
        
        shutil.copytree(file, ROOTDIR + file)
            
if __name__ == "__main__":    
    init()
    check_python()
    check_CUDA()
    INSTALLDIR = get_installdir()
    ROOTDIR = create_installdir(INSTALLDIR)
    install(ROOTDIR)