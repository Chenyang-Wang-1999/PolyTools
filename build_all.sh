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
cp python/py_invariant_manifold.py build/py_invariant_manifold/__init__.py
cp python/*.so build/py_invariant_manifold/
cp test/*.py build
