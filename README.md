# Tycho2 Version 0.2

A mini-app for neutral-particle, discrete-ordinates (SN), transport on parallel-decomposed meshes of tetrahedra.

## Background and History

Tycho2 is based on a code called Tycho, written by Shawn Pautz (Sandia National Laboratory) sometime around the year 2000,
when he was at Los Alamos National Laboratory. Though Tycho2 wasn't a mini-app (that term hadn't been invented yet), and 
didn't borrow any code base from the original implementation, Shawn was using Tycho as a prototyping code to test ideas for
parallel, unstructrued mesh, SN sweeps. Remember, those were early days. We named our mini-app Tycho2 in an homage to Shawn's 
original code (and because Tycho2 works only on tets, like Tycho).

Tycho2 was originally written by Kris Garrett while a postdoc at Los Alamos National Laboratory (see the license) in the
CCS-2 group, working together with Jim Warsa and Jae Chang, members of the Capsaicin deterministic SN transport code project 
at LANL. The goal of the mini-app is to provide a simple platform to explore the strange new world of GPUs and heterogeneous
computer architectures. The idea being that, like the SNAP, UMT, and Kripke mini-apps for structured meshes*, potential approaches 
and implementations for unstructured mesh transport could be investigated and perhaps transferred to the more complicated Capsaicin
project if they showed promise. Makes sense, right? Anyway, Kris has left LANL and now it's up to Jim Warsa and any other contributors
to take up where he left off. Recent work by Dan Ibanez of Sandia National Laboratory, with Kris, on using Kokkos in Tycho2 is
gratefully acknowledged.

*Note that, to the best of our knowledge, while UMT is targeted to unstructred meshes, the meshes are actually logially 
orthgonal, without the additional complications and computational overhead associated with SN sweeps on meshes with 
arbitrary orientations and connectivtiies.

## Warning
This is still being heavily developed, so ...
USE AT YOUR OWN RISK!!!


## License
Copyright (c) 2016, Los Alamos National Security, LLC  
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:  
1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.  
2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.  
3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Build Notes
Build notes can be found at the following wiki page:
https://github.com/losalamos/tycho2/wiki/Build-Notes


## Los Alamos LACC Number
LA-CC-16-049


## Contributors
- Kris Garrett (LANL [formerly], original author)
- Neelam Patel (Affiliation unknown)
- Kevin Procopio (ORNL, contributor)
- Jim Warsa (LANL, ideas, criticisms, advisement)
- David Dixon (LANL, contributor)
- Dan Ibanez (SNL, contributor, Kokkos expert)
- Jae Chang (LANL, honorable team leader)
