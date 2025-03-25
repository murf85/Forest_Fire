# Forest_Fire

Cell-DEVS Model of a forest fire based on the Corsica Model

Introduction
This repository contains a Cell-DEVS model of an Aircraft Transponder.  In addition to the top model, there are two coupled models representing a radar and an aircraft. 
The radar sends interrogation pulses to the aircraft, and the aircraft replies with Mode A code or Mode C altitude.   

Dependencies
This project assumes that you have Cadmium installed in a location accessible by the environment variable $CADMIUM. 

Build
To build this project, run:

source build_sim.sh

NOTE: Every time you run build_sim.sh, the contents of build/ and bin/ will be replaced.

Execute
In order to run the model a .json configuration file is required. The examples discussed in the report can be found in the config folder.  

To run the models in this project, replace the <filename> with name of .json file containing the simulation parameters and input:

./bin/forestfire config/<filename>.json

The program will output a log file. Log files are best explored using the Cell-DEVS viewer found here:

https://devssim.carleton.ca/cell-devs-viewer/

Both the config and log files are required for the viewer to operate.  


