# Alenka-Signal

This repo has been merged with the parent project [Alenka](https://github.com/machta/Alenka).

### Requirments
* git
* cmake 3.1 (for Ubuntu 14 download the latest version)
* g++ or Microsoft Visual C++ 2015 compiler
* OpenCL 1.1

On Debian-like systems you can use: `sudo apt install git cmake-gui build-essential`

MSVC compiler can be acquired by installing Visual C++ Build Tools 2015. Choose Custom
installation, and uncheck all options but "Windows 8.1 SDK". If you already have Visual
Studio, you don't need to install this.

If your device doesn't support OpenCL (e.g. when running a Linux guest in VirtualBox),
use AMD APP SDK for a CPU implementation of OpenCL.

### Build instructions
1. Clone this repo
2. Install libraries using download-libraries.sh
3. Make new build directory and change into it
4. Use cmake-gui to generate build environment
   1. Click Configure and choose a compiler
   2. Change CMAKE_BUILD_TYPE (not needed on Windows)
   3. Change CMAKE_INSTALL_PREFIX
   4. For a 32-bit build make sure BUILD64 is unset
   5. Click Generate
5. Build the library

Here is an example of the setup using git-bash (or regular bash):
``` bash
./download-libraries.sh
mkdir build-Release-64 && cd build-Release-64
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install-dir ..
# On Windows: cmake -D CMAKE_INSTALL_PREFIX=install-dir -G "Visual Studio 14 2015 Win64" ..
cmake --build . --config Release --target install-alenka-signal
```

