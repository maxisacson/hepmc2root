# hepmc2root
Convert from the HepMC2 format to a flattned ROOT tree.

## Build instructions
Dependencies: CMake, HepMC2, ROOT6
```
$ mkdir build
$ cd build
$ cmake ../ <CMake options>
$ make
```

CMake options:
- `-DHEPMC_PATH=<path>` Path to the HepMC2 installation. Default: `/usr/local`
- `-DROOT_CONFIG=<path>` Path to the `root-config` executable if not in $PATH. Default: `root-config`
- `-DCMAKE_BUILD_TYPE=<Debug|Release>` Build type. Default: `Debug`

## Usage
```
$ ./hepmc2root input.hepmc
```

For more info see
```
$ ./hepmc2root -h
```
