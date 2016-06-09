# Tycho 2 Version 0.1

Do not use yet.
Still creating the repository.

## Typical Usage

In the `util` directory, there are several mesh files already created.
You will need to make a utility to partition a mesh into parallel meshes.
You can use either `SerialToParallelMesh` or `PartitionColumns`.
This will create a `.pmesh` file for the sweep program to use.

## Building the main program

To build the main program, copy `make.inc.example` to `make.inc` and fill in your local parameters.
Then type `make` to build the program, sweep.x.
