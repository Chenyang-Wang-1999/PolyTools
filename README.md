# PolyTools

A C++ package for efficient polynomial calculations.

```
Note: current version is based on a pure C linked list, which is unefficient. Besides, I used this project for both the nonlinear dynamics and the non-Hermitian bands, so there are different python packages other than "poly_tools". I will refactor the whole project if I have time. 

But now, I am struggling for a job and could not spare any time. So, please pardon me for leaving such a mess in this project.
```

Prerequisite: You need pybind11 installed compatible to your python version.

Compile:
``` bash
cd ./pybind
make _poly_tools_cc.<your-python-lib-suffix.so>
```

Then, you can copy the folder `./python/poly_tools` to your working directory and use it.

Apologize again for such a s**t.
