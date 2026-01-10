## Molecular Dynamics Simulation Project  
                    
This project implements a Molecular Dynamics (MD) simulation in Fortran 90, utilizing the Lennard-Jones potential and the Verlet integration algorithm to study atomic trajectories.

### Directory Structure

The repository is organized as follows:
* **src**: Contains the source code file `dynamics.f90`.
* **tests**: Contains an example of the needed imput files and of the expected output.
    * param.inp: Input file for simulation parameters.
    * inp.txt: Input file containing the initial coordinates and masses of the atoms.
    * trajectory.xyz: Output file that is expected when using those input files.* inp.txt:
* **README.md**: Describes the project, its physical background, and structure.
* **INSTALL.md**: Contains instructions on how to compile and execute the code.
* **AUTHORS**: Lists the contributors to this project.
* **LICENSE**: Specifies the licensing terms for the project.

### Program Overview: dynamics.f90

#### Input files
The input files can have the same shape of the test ones.
* **param.inp**:
    * epsilon must be espressed in kJ/mol.
    * sigma must be expressed in Å.
    * deltat time step expressed in ps.
    * Nstep numer of the time steps.
    * symbol is the name of the atom.
* **inp.txt**:
    * Natoms as first line.
    * one line each atom containing the cartesian coordinates expressed in Å and mass espressed in g/mol.

#### Main, subroutines and functions
* The main program allocates variables and arrays after reding the number of atoms.
* It reads the Lennard-Jones parameters sigma and epsilon and the one needed to the simulation from `param.inp`.
* The `read_molecule` subroutine stores the masses and initial coordinates of each atom.
* The `compute_distances` subroutine calculates the distance between every pair of atoms, storing them in a 2D array.
* The `compute_acc` subroutine determines the acceleration of each atom based on the Lennard-Jones potential.
* Energy Functions: Potential (V), Kinetic (T), and Total Energy (E)
* `Verlet` subroutine

#### The Verlet Algorithm
The `Verlet` subroutine implements the integration process to update positions and velocities at each time step.
1. Starting velocities are set to zero.
2. At each step, the subroutine calls the compute_distances and compute_acc subroutines to update distances,accelerations and
   velocities.
3. At each step, the energy functions are called to uptade the energy components values.
The simulation ends when the total number of iterations reaches `Nstep`.
The parameter needed to the Verlet algorithm can be changed from the param.inp file, whose structure was previously defined.

### Output

The program generates a `trajectory.xyz` file, writing every **10 steps**:
1. The total number of atoms.
2. A comment line Containing Kinetic, Potential, and Total Energy values to verify energy conservation.
3. A line for each atom including the chemical symbol followed by its X, Y, and Z coordinates in Å.

The animation of the trajectory can be visualized with molden by using:

```bash
molden dynamics.f90 
```
Once inside molden you can press `movie` to see the trajectory and `Geom. conv.` to see if the energy is conserved.
