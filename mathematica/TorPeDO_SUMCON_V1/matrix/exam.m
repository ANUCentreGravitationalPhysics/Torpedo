matKpvc = {{-505.03154719946906 - 25.237197706428788*I}}
 
matKvvc = {{505.03154719946906 + 25.237197706428788*I}}
 
matMvv = {{320}}
 
matMpv = {{0}}
 
matGpv = {{0}}
 
matGvv = {{0}}
 
allvars = {yPL}
 
allvars /: allvars::usage = "allvars is part of the specification of the \
model. It should be a list of the position and angle variables for all \
elements of the system that have mass/MOI associated with them."
 
allparams = {yg}
 
allparams /: allparams::usage = "allparams is part of the specification of \
the model. It should be a list of the positions and angles describing the \
state of the ground."
