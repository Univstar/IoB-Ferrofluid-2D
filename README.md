# IoB-Ferrofluid (2D)

This repository maintains the code of the SIGGRAPH 2024 paper *An Induce-on-Boundary Magnetostatic Solver for Grid-Based Ferrofluids*.

Here we show a 2D example of ferrofuild simulation, which is not accelerate by FMM because the number of points is relatively low.
This program can be performed on Windows, Linux, and MacOS, with one of the following compilers installed:
- MSVC 19.30+
- GCC 11+
- Clang 13+
- Apple Clang 13+

To play with the code, please first install [xmake](https://xmake.io/) as instructed in https://xmake.io/#/getting_started.
Then, run the following commands successively in this directory for getting started:
```shell
xmake
xmake r demo -t box -r 500 -e 801 -s 256
```
The exported images are subsequently generated in `build/[Platform]/[Arch]/release/output`.

We acknowledge [the work](https://jcgt.org/published/0011/02/02/) of Tetsuya Takahashi and Christopher Batty for [MC-style-vol-eval](https://github.com/tetsuya-takahashi/MC-style-vol-eval).
