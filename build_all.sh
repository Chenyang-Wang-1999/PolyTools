#!/bin/bash

# if ! [ -e build ] ; then
#     mkdir build
# fi

rm -r build
rm obj/*
rm pybind/*.so

mkdir build
mkdir build/py_invariant_manifold

make all
cd pybind
make all
cd .. 

cp pybind/*.so python/py_invariant_manifold
cp -r python/* build