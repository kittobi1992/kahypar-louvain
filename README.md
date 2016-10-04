KaHyPar - Karlsruhe Hypergraph Partitioning
=========
Travis-CI Status: [![Travis-CI Status] (https://travis-ci.com/SebastianSchlag/kahypar.svg?token=ZcLRsjUs4Yprny1FyfPy&branch=master)](https://travis-ci.com/SebastianSchlag/kahypar.svg?token=ZcLRsjUs4Yprny1FyfPy&branch=master)
Appveyor Status: [![Appveyor Status](https://ci.appveyor.com/api/projects/status/5kxoc64byyhae8ue/branch/master?svg=true)](https://ci.appveyor.com/api/projects/status/5kxoc64byyhae8ue/branch/master?svg=true)
Coverity Status: [![Coverity Status](https://scan.coverity.com/projects/10262/badge.svg)](https://scan.coverity.com/projects/10262/badge.svg)

Requirements:
-----------
The Karlsruhe Hypergraph Partitioning Framework requires:

 - A 64-bit operating system. Both Linux and Mac OS X are currently supported.
 - A modern, C++11 ready compiler such as `g++` version 4.9 or higher or `clang` version 3.2 or higher.
 - The [cmake][cmake] build system.
 - The [Boost.Program_options][Boost.Program_options] library.


Building KaHyPar:
-----------

1. Clone the repository including submodules: `git clone --recursive git@github.com:SebastianSchlag/kahypar.git`
2. Create a build directory: `mkdir build && cd build`
3. Run cmake: `cmake .. -DCMAKE_BUILD_TYPE=RELEASE`
4. Run make: `make`

Testing and Profiling:
-----------

Tests are automatically executed while project is built. Additionally a `test` target is provided.
End-to-end integration tests can be started with: `make integration_tests`. Profiling can be enabled via cmake flag: `-DENABLE_PROFILE=ON`.  

Running KaHyPar:
-----------

The binary is located at: `build/kahypar/application/`.

KaHyPar has several configuration parameters. For a list of all possible parameters please run: `./KaHyPar --help`
    
Currently we provide two different presets that correspond to the configuration used in the 
[ALENEX'16](http://epubs.siam.org/doi/abs/10.1137/1.9781611974317.5) submission and the [ALENEX'17]() submission.

To start KaHyPar in recursive bisection mode optimizing the cut-net objective run:

    ./KaHyPar -h <path-to-hgr> -k <# blocks> -e <imbalance (e.g. 0.03)> -o cut -m recursive -p ../../../config/cut_rb_alenex16.ini
    
To start KaHyPar in direct k-way mode optimizing the (connectivity - 1) objective run:   
  
    ./KaHyPar -h <path-to-hgr> -k <# blocks> -e <imbalance (e.g. 0.03)> -o km1 -m direct -p ../../../config/km1_direct_kway_alenex17.ini
    
All preset parameters can be overwritten by using the corresponding command line options.


Bug Reports:
-----------

We encourage you to report any problems with KaHyPar via the [github issue tracking system](https://github.com/SebastianSchlag/kahypar/issues) of the project.


Licensing:
---------

KaHyPar is free software provided under the GNU General Public License (GPLv3). 
For more information see the [COPYING file][CF].

We distribute this framework freely to foster the use and development of hypergraph partitioning tools. If you use KaHyPar in an academic setting please cite the following paper:
    
    @inproceedings{shhmss2016alenex,
     author    = {Sebastian Schlag and
                  Vitali Henne and
                  Tobias Heuer and
                  Henning Meyerhenke and
                  Peter Sanders and
                  Christian Schulz},
     title     = {\emph{k}-way Hypergraph Partitioning via \emph{n}-Level Recursive
                  Bisection},
     booktitle = {18th Workshop on Algorithm Engineering and Experiments, (ALENEX 2016)},
     pages     = {53--67},
     year      = {2016},
    }

A preliminary version is available [here on arxiv][ALENEX16PAPER].

Contributing:
------------
If you are interested in contributing to the KaHyPar framework
feel free to contact any of the authors or create an issue on the
[issue tracking system](https://github.com/SebastianSchlag/kahypar/issues).


[cmake]: http://www.cmake.org/ "CMake tool"
[Boost.Program_options]: http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html
[ALENEX16PAPER]: https://arxiv.org/abs/1511.03137
[CF]: https://github.com/SebastianSchlag/kahypar/blob/master/COPYING "Licence"
