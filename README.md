# Design of steel-timber composite beams #

This is simple python code written to replace the manual computations for 
the design of the slim-floor steel-timber composite beams (Nordic system). 
The composite beam includes a WQ-beam with two CLT slab resting on its 
lower flanges.
The composite action is realized by connecting the slabs to the WQ-beam 
by steel plates. 

The analytical background is developed in the manuscript 
"Effective width of slim-floor steel-timber composite beams" by 
Garifullin, Mela, Pajunen, Aspila and Heinisuo submitted for publication in
the Journal of Constructional Steel Research.

Any publication of this code without the permission of the author 
is not allowed.

### How to use this code? ###
* In wq_clt_input.py, set the input properties 
* In wq_clt_main.py, set the coordinates of the points at which 
the desired parameters to be calculated and assign them to the requested 
functions
* Run wq_clt_main.py
* The requested results are printed to the console

### Realized functions ###
* composite_beam.vertical_displacement(x): vertical displacement of the WQ-beam
* composite_beam.stresses_wq(x): axial stresses in the top and bottom 
flanges of the WQ-beam
* composite_beam.stresses_clt(x): stresses in the CLT slab in the direction 
of the WQ-beam
* composite_beam.screw_force(x): force in the screw of the connector
* composite_beam.vibration(): fundamental structural frequency of the 
composite beam

### Contacts ###
Marsel Garifullin: marsel.garifullin@tuni.fi 

Kristo Mela: kristo.mela@tuni.fi

