# ECH 267 Advance Process Control: Final Project 

# Introduction: 

The goal of this project is to design, develop, and implement a MPC controller and path planner which is capable of driving a simple two link robot from one pose to another. Ideally, this would include avoiding ** static ** obstacles which obstruct the controller from simplest and direct path, as well as obeying constraints on the permitted range of motion which the system is subject to. 


# Directions: 

To run the scripts in the following `2_Code` folder, please follow the instructions below. 

1. Install Matlab 
2. Download and install Casadi from the instructions below. 
3. Clone this repository to a folder on your local PC. 
4. Add this folder to the Matlab path. 
5. Run scripts


# System Details: 
The following are the specifications of the system which I have will be using to implement this project. 

* OS: Windows 10 
* MATLAB Version: R2020a  
* CasADi Version: v3.5.5


## Optimization Framework: 

CasADi, is an opensource numerical/symbolic optimization and automatic differentiation package which enables a streamlined matlab or python implementation of complex optimization problems. It has the desired effect of abstracting the complicated optimization algorithms into an easy to use syntax which in turn is used to define the optimization problem symbollically. This reduces the size of boiler-plate code as well as enabling this project to focus on the MPC and not the optimization details. 

## Installing CasADi:

The language which this project uses is MATLAB as there are variety of well documented implementations of the control, optimization, and visualization elements which are not nearly as straightforward in python, in the authors opinion. 

To install CasADi ... 
1. Click the following link: https://web.casadi.org/get/
2. Download the package which fits your system 
3. Unzip the downloaded folder 
4. Add the unzipped folder to the MATLAB Path. 
5. Include `import casadi.*` if not already written in the script





