Fluid-Structure Interaction (FSI) solver using the Fictitious Domain Method.

The finite element formulation is published in the paper in [CMAME](https://www.sciencedirect.com/science/article/abs/pii/S0045782515004296) journal.
The staggered scheme is the one in this [paper](https://scholar.google.co.uk/citations?view_op=view_citation&hl=en&user=zNk33SUAAAAJ&citation_for_view=zNk33SUAAAAJ:qxL8FJ1GzNcC).

* Programming languages: **C++** and **Fortran**
* C++ standard: **C++14** (or above)
* Required third-party libraries:
  1. CMake
  2. Blas
  3. Lapack
  4. Boost
  5. MPI (OpenMPI, MPICH or Intel MPI. Your choice!)
  6. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  7. [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
  8. [PETSc](https://www.mcs.anl.gov/petsc/)
  9. [VTK](https://vtk.org/)
  10. Metis and ParMetis
  11. SuperLU


## Compilation and building
1. Clone the repository or download the zip file and extract its contents.
2. Go to the directory of the repository in a terminal.
3. Create **build** and **bin** directories.
    * `mkdir build`
    * `mkdir bin`
4. Modify the CMake file accordingly.
    * Copy `CMakeLists-chennalaptop.txt` to a new file for your machine, say `CMakeLists-local.txt`.
    * Change the paths to the compilers, Eigen, CGAL, PETSc and VTK libraries.
    * Change the path in the `install` function.
    * Create a symbolic link to the local CMake file.
      * `ln -sf CMakeLists-local.txt CMakeLists.txt`
5. Enter the `build` directory.
    * `cd build`
6. Configure using the CMake file
    * `cmake ..`
7. Compile, build and install the executable `fsifdm`. This step will also copy the exe to the `bin` folder.
    * `make install`

## Execution
* Simulations are usually run from the `bin` folder.
* For serial run: `./fsifdm  <pathtodirectory>  <inpfilename>` (White space between each entry)
* Example:
  * `./fsifdm  ../project/sampleinputs  IGallopingFDM`

## Understanding and using output(s)
For each simulation the output contains two items:

1. Files for visualisation in ParaView. Generated in the same directory of the input file.
    * One **.vtu** file for each time step. Contains velocity, pressure and vorticity for the flow field.
    * One **.vtu** file for each solid, if the solids are allowed to move. Only done for flexible solids.

2. A file with the letter prefixed `T` to the name of the input file containing the data for forces and displacements. (This will be created in the same directory as the input file.)