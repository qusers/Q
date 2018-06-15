---
project: Q-FEP-MD
summary: Developer docs for Q
author: John Marelius, Johan Aqvist, and others.
author_description: Research group at Uppsala
project_dir: ../../../src 
output_dir: ./devdocs
media_dir: ./media
project_github: https://github.com/qusers/qsource
project_website: http://qdyn.no-ip.org
source: false
graph: true 
macro: TEST LOGIC=.true.
github: http://www.icm.uu.se/cbbi/aqvist-lab 
email: johan.aqvist@icm.uu.se  
predocmark: >  
predocmark_alt: <  
docmark: !  
coloured_edges: true  
sort: alpha  
display: public protected private  
---

--------------------


Introduction {#introduction}  
============

This  is  the  developers  manual  for **Q**,  the  free  energy  from
molecular dynamics program.  

**Q**  is  entirely  written   in  modern  **FORTRAN**.   Its  initial
development was  in **FORTRAN**  77, but  the incorporation  of object
oriented paradigms into  modern **FORTRAN** have seen  the code evolve
to take  advantage of  the evolution of  the language.   Currently the
program  incorporates elements  of the  **FORTRAN** 2008  standard and
it's developing the use of the OpenMP 4.0 standard (*fully implemented
into gnu Fortran version 5.1*)
for in-node parallelization.  


License  
=======

For a  non-commerical license agreement please  fill-out the following
form and send it
to Johan Aqvist.  

[LICENSE](http://www.icm.uu.se/digitalAssets/211/211337_3q_license.pdf)  

In return you will be given access to the code.  

