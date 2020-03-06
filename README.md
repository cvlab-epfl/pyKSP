### Python wrapper for “Multiple Object Tracker Using K-Shortest Paths”
Wrapper written by Pierre Baqué from the original "Multiple Object Tracker Using K-Shortest Paths" code.

### Usage
If you want to modify the source code or need recompile the wrapper, enter the command:
```
cd pyKSP
python setup.py build_ext --inplace
```

Or install to your site-packages to make the module available everywhere:
```
cd pyKSP
python setup.py install
```

The python script in "pyKSP-Example.py", provides an example on how to use pyKSP from Python.

### Dependencies

Python libraries:
```
pip install cython numpy matplotlib
```

[Boost Graph Library](https://www.boost.org/doc/libs/1_72_0/libs/graph/doc/index.html) to build the KSP binary:
```
# on ubuntu
sudo apt install libboost-graph-dev
```

### References
For more information about the KSP algorithm, please check the
following article:

Jerome Berclaz, Francois Fleuret, Engin Turetken and Pascal Fua,
"Multiple Object Tracking using K-Shortest Paths Optimization", IEEE
Transactions on Pattern Analysis and Machine Intelligence (TPAMI),
2011.
