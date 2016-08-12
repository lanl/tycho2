# Tycho 2 Version 0.1


## Warning
This is still being heavily developed.
Use at your own risk!!!


## License
Copyright (c) 2016, Los Alamos National Security, LLC  
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:  
1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.  
2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.  
3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Quick start guide

The first step is to create a parallel mesh.
To do this:
- Go to the `util` folder
- Copy `make.inc.example` to `make.inc`
- In `make.inc`, write in your C++ compiler with the option for the C++11 standard set.
For example with GCC, this would be: g++ -std=c++11.
- At the command prompt, type `make PartitionColumns`
- Then type `./PartitionColumns.x` to see how to use the utility.
- Type `./PartitionColumns.x 1 1 cube-1374.smesh cube-1374.pmesh` to create a mesh that is only partitioned into 1 partition.
- Go back to the previous directory `cd ..`

Now that you have created a parallel mesh, it is time to compile the main program.
- Copy make.inc.example to make.inc
- Copy input.deck.example to input.deck
- In `make.inc`
  - Set ASSERT_ON = 0
  - Set MPICC = <your MPI compiler with C++11 standard and OpenMP> (e.g. mpicxx -std=c++11 -fopenmp)
- Type `make`.  This will build `sweep.x`
- Then run `./sweep.x util/cube-1374.pmesh input.deck`


## Los Alamos LACC Number
LA-CC-16-049


## Contributors
- Kris Garrett (Original Author)
- Neelam Patel
- Kevin Procopio
