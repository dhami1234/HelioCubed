# SWQU User's Meeting Tutorial

In this part of the tutorial, we will use the cubed sphere grid to create an SW solution using a time-independent BC generated by SWiG. The BC has only one frame and it is rotated with the solar rotation rate at the inner boundary.

The BC file has been placed in the `/BC_files` directory and it is also added to the inputs file we will use when running the code.

## Steps for Running the Code

1. **Navigate to the proto directory and set the branch to `cubedSphereDev` using:**
   ```bash
   git branch cubedSphereDev

2. **Navigate to HelioCubed directory/exec.**
   ```bash
   cd exec

3. **Edit the GNUMakefile with the correct paths to proto, openblas, and hdf5 and save the file.**

4. **Compile HelioCubed using**
   ```bash
   make -j8

5. **The code can now be run using**
    ```bash
    mpirun -n 6 CubedSphere.exe inputs

6. **The inputs file can be edited to set desired parameters of the run.**


