# ECH 267 Advance Process Control: Final Project 

#Introduction: 

The goal of this project is to design, develope, and implement a MPC based Trajectory-Generator for a Robot arm. The desired outcome of this project is to use MPC on a nonlinear discrete dynamics model of the robot and generate a trajectory for it to follow using the principle of optimal contol as applied via MPC. The secondary goal of this project is to further develope a controller capable of following the 'optimal' trajectory it is being given. Time permitting it would get to actually develop a sensor model for this robot and include obstacles so tha the trajectory might need to change "in close to real-time" as soon as the sensor detects an obsticle. 

# Directions: 



# System Details: 
The following are the specifications of the system which I have will be using to implement this project. 

* OS: Windows 10 
* MATLAB Version: R2020a  
* CasADi Version: v3.5.5


## Optimization Framework: 

CasADi, is an opensource numerical/symbolic optimization and automatic differentiation package which enables a streamlined matlab or python implementation of complex optimization problems. It has the desired effect of abstracting the complicated optimization algorithms into an easy to use syntax which in turn is used to define the optimization problem symbollically. This reduces the size of boiler-plate code as well as enabling this project to focus on the MPC and not the optimization details. 

## Installing CasADi:

The tentative language which this project will most likely use is MATLAB as there are variety of well documented implementations of the control, optimization, and visualization elements which are not nearly as straightforward in python. 

To install CasADi ... 
1. Click the following link: https://web.casadi.org/get/
2. Download the package which fits your system 
3. Unzip the downloaded folder 
4. Add the unzipped folder to the MATLAB Path. 
5. Include `import casadi.*` if not already written in the script





