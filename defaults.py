MD = {  
'shake_solvent'     : False
} 

cutoffs = {
'solute_solvent'    : 10    ,
'solute_solute'     : 10    ,
'solvent_solvent'   : 10    ,
'q_atom'            : 99    ,
'lrf'               : 99    ,
} 

sphere = {
'shell_force'       : 10.0  ,
'shell_radius'      : 25    ,
} 

solvent = {
'radial_force'      : 60.0  ,
'polarisation'      : True  ,
'polarisation_force': 20.0  ,
} 

intervals = {
'output'            : 5     ,
'trajectory'        : 100   ,
'non_bond'          : 25    ,
'energy'            : None  ,
} 
