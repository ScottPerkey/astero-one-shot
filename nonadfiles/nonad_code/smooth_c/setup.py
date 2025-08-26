#!/usr/bin/env python

# this is a generic installation file for a python code, taken from Daniel Foreman-Mackey.

import re
import os
import sys
from setuptools import setup
from Cython.Build import cythonize
from distutils.core import Extension
from setuptools import setup, find_packages
# for cython:
# compile with:
# python setup.py build_ext --inplace
name = 'smooth_c'
package_names = [name]

try:
    version = re.findall(r"__version__ = \"(.*?)\"", open(module+".py").read())[0]

except:
    version = '0.0'

setup(
    name=name,
    version=version,
    author="Joel C. Zinn",
    author_email="zinn.44@osu.edu",
#    url="https://github.com/dfm/triangle.py",
    package_dir={'smooth_c': ''},
    packages=package_names,
#    package_data={'/home/ritchey/zinn.44/kochanek/test_prh/':["qso_periodic.cfg", "rrlyr_periodic.cfg", "lpv_periodic.cfg", "qso_flicker.cfg", "rrlyr_flicker.cfg", "lpv_flicker.cfg", "mcmc.cfg"]},
#    data_files = [(name, ["mcmc.cfg"])],
#    package_data={'': package_data_files},
#    include_package_data = True,
    # ext_modules=[Extension('band_construct_c', ['band_construct_c.so']), Extension('band_matrix', ['band_x_vec.so'])]
# don't want to assume that cython is installed in the install location. so will just build a .so file and attach it to the distribution.    
#    ext_modules = cythonize(cython_modules)
#    description="",
#    long_description=open("README").read(),
#    package_data={"": ["LICENSE"]},
#    include_package_data=True,
#    classifiers=[
#        "Development Status :: 5 - Production/Stable",
#        "License :: OSI Approved :: BSD License",
#        "Intended Audience :: Developers",
#        "Intended Audience :: Science/Research",
#        "Operating System :: OS Independent",
#        "Programming Language :: Python",
#    ],
)


