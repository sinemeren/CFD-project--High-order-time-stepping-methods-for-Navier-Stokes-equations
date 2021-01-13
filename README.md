# CFD Lab Group Project #

This repository was created for the *Computational Fluid Dynamics Lab - Summer 2020* course group project.

Authors: Berkay Sezen, Sultan Sinem Eren and CFD Lab code skeleton contributors

# Project Description #

In this project, different fixed time-stepping methods were implemented to solve various problems that perform *Incompressible Navier-Stokes* equations. The purpose of the implementation is to investigate the behavior of different numerical methods in both stationary and transient cases with varying time-steps. The measurement criterion is based on simulation time, number of iterations number and solution accuracy. Implemented methods and scenarios are shown below. 

### Problems ###
  * Lid-Driven Cavity
  * Plane Shear Flow
  * Flow Over a Step
  * The Karman Vortex Street
  
### Time-Stepping Methods ###
  * Explicit Euler
  * Heun (Runge-Kutta 2nd order)
  * Runge-Kutta 4th order

# Code Description #

This repository contains

- Source code for the simulations
- Source folder with required source files (*.dat* and *.pgm* files)
- Additional tester and documentation folders and codes

In order code to be run, at first, it should be compiled. Steps for compiling are as follows:

1. Create a build directory: `mkdir build`
2. Get inside it: `cd build`
3. Configure and generate the build system: `cmake ..` (Note the two dots, this means that the `CmakeLists.txt` File is in the folder above)
4. Build your code: `make` (build the executable)

## Running a Simulation ##

There are multiple options for both problems and methods. Therefore, please use the following format: `./sim [problem] [method]`

The available problems are:
- LidDrivenCavity
- FlowOverStep
- PlaneShearFlow
- KarmanVortexStreet

The available methods are:
- ExplicitEuler
- Heun
- RungeKutta

An example simulation command is:
```shell
 ./sim LidDrivenCavity ExplicitEuler
```

Sample output at command window is : `t = 2.060000 ,dt = 0.005000, Res = 0.000998,iterations=18`

After running a simulation, respective folders and files are created under `Result/` folder according to used *problem* and *method*. Created *.vtk* files can be used in *Paraview* software for visualization purposes. Also, three different *.txt* files are created to save `number of iterations`, `residuals` and `corresponding time values`. It is possible to use these *.txt* files for comparison and plotting purposes.

## General Remarks ##

It is also possible to run a simulation with different configurations and scenarios.
- Step-size can be adjusted inside main.cpp between lines 221-227 (requires re-compiling). It also is possible to use both fixed-time steps and adaptive time steps.
- Configuration files for the problems can be found under `Source/` folder. Implementation is independent of geometry so that by changing *.pgm* files inside `Source/pgm/` different scenarios can be observed. Also, the parameters for the simulations (*.dat* files) can be found under `Source/dat` so that simulations with different parameters are possible.

# Additional Information #

## Software Requirements

* VTK 7 or higher
* GCC 9 (optional) 

### GCC version

You can get you current version of GCC by running:

```shell
g++ -v
```

### Defining your GCC version

If you have GCC 9 or newer, you can set in the `CMakeLists.txt` file:

```cmake
set(gpp9 True)
```

If you have a version lower than 9, then you don't have to modify the `CMakeLists.txt` file.

This will affect how we are using the C++ filesystem library, which is available already in GCC 7 as an experimental feature.

### Setup of VTK and GCC 9 (Ubuntu **20.04**)

```
apt-get update &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev
```

### Setup of VTK and GCC 9 (Ubuntu **18.04**)

If you want, you can upgrade your compiler version to have access to more recent C++ features.
This is, however, optional.

```
apt-get update &&
apt-get install -y software-properties-common &&
add-apt-repository -y ppa:ubuntu-toolchain-r/test &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev gcc-9 g++-9
apt-get install -y gcc-9 g++-9
```

### Troubleshooting: VTK not found

You might run into a problem where the VTK library is not found. To fix this, you can try the following steps:

1. Find the installation path of your VTK library 
2. Define this path as an environment variable, as e.g. `export VTK_DIR=".../lib/cmake/vtk-8.2"`
3. Start in a clean build folder
4. Run `cmake ..` again

### Set a different GCC version

If you have multiple compiler versions installed you can set the GCC version which should be used by `cmake` like this:

```shell
export CXX=`which g++-7`
```

Make sure to use a backtick (\`) to get the `which` command executed. Afterwards, you can run `cmake ..`.

## Testing 

The folder `tests` contains example files for simple unit testing with the [Catch2](https://github.com/catchorg/Catch2) unit testing framework. The tests can be compiled using `cmake` as well.

After building, run:

```
ctest --verbose
```

With the `verbose` option you can get more details for failing tests.

## Documentation 

The generate the documentation, run:

```
pip3 install --user mkdocs mkdocs-material
```

and then serve it to be able to access it through your browser:

```
mkdocs serve
```

