# TMVA_BDT_FRAMEWORK

# step 1 : Make samples 
Prepare your input samples in the format of particle-by-particle.
For instance, if you have 10 particle candidate in a single collision event, your 'input tree' should be designed to have 10 events.
Make sure the format is identical for signal and background

-- You might need to assign a 'weight' for signal and background. Even if you want to weight only signal events, put the same branch also into backgrounds with a number equal to 1.


# step 2 : BDT Training 
The main training part - run BDTClassifier.C
You need to define 
1) Input files for signal and background / output file 
2) Training variables 
3) Spectators (not training variables)
4) Pre-selections 
5) number of trees, number of events to train, Max Depth, Min Node size, Boost Type, Boost Beta, Seperation type, Pruning, RandomTrees, nCuts, etc. 
-- mostly just start with the default setting 

# step 3 : BDT application
This is for applying the trained BDTs into the trees you want to apply - run BDTClassifierApplication.C
You need to define 
1) BDT xml file you want to use (outcome of step 2)
2) Input files you want to apply / Output file
