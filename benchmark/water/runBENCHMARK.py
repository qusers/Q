import os
import shutil

class Create_Environment(object):
    """
        Creates the workdirectory environment.
    """
    def __init__(self,top,wd):
        if not os.path.exists(wd.split('/')[0]):
            os.mkdir(wd.split('/')[0])
            
        else:
            shutil.rmtree(wd.split('/')[0])
            os.mkdir(wd.split('/')[0])
        
        if not os.path.exists(wd):
            os.mkdir(wd)
            
        else:
            shutil.rmtree(wd)
            os.mkdir(wd)
            
        # create the output folder for qdyn
        os.mkdir(wd + '/output')

class Prepare_Q_topology(object):
    """
        Prepares a Q topology using qprep (Q6.0.7 Fortran code)
        
    """    
    def __init__(self):
        self.qprep_inp = """ """"


class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qdyn:
               {'top'       :   top,
                'fep'       :   fep,
                'md'        :   md,
                're'        :   re,
                'wd'        :   wd,
                'verbose'   :   verbose
                'clean'   :   clean
               }
        """

        # INIT
        #Create_Environment(top = self.environment['top'],
        #                   wd  = self.environment['wd'],
        #                )