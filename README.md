# ZoomPM - Super resolution with Perona Malik Diffusion

Version 0.9, 07/04/2024

## Authors
- Abdelmounim Belahmidi <abelahmidi@lerity-alcen.com>
- Pascal Monasse <pascal.monasse@enpc.fr>

## Licensing
See LICENSE file

## Build
The build procedure is handled by CMake (https://cmake.org/).
Required library: [libPNG](http://libpng.org/pub/png/libpng.html)

Under Linux, the following commands are enough:
```
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release /path/to/source
$ cmake --build .
```

## Usage
```
./zoomPM [options] imgIn.png imgOut.png
Options:
-z ARG Integer zoom factor (2)
-c, --color Handle color image
-n, --iter=ARG Nb of iterations (0=auto) (0)
```

The zoom factor is an integer.

By default, color images are converted to grayscale; use -c option to force handling the three channels.

The number of evolution iterations is automatically computed (and displayed) by default. It is not useful to increase this number, but you can divide it by 2 without dramatic loss of quality.

## Example
Try with the images of folder data, for instance:
```
$ ./zoomPM -z 4 -c ../data/zebra.png zoom.png
```
Compare image zoom.png with ../data/zebra_x4.png.

## Files
- CMakeLists.txt, LICENSE, README.md
- folder data/
- cmdLine.h, io_png.{c,h}, xmtime.h: support code
- zoomPM.cpp (IPOL reviewed)
