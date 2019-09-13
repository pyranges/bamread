import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

__version__ = open("bamread/version.py").readline().split(" = ")[1].replace(
    '"', '').strip()
macros = []

install_requires = ["scipy", "numpy", "natsort", "cython", "pysam", "pandas"]

compile_options = [
    "-Ofast", "-Wall"
    # , "-stdlib=libc++"
]  #, "-frename-registers", "-funroll-loops"] # , "-lgzstream", "-lz"

dir_path = os.path.dirname(os.path.realpath(__file__))
lib_dirs = []
include_dirs = [dir_path + "/bamread/src", dir_path]

from subprocess import check_output

# conda_path = check_output("which conda", shell=True).decode().strip()
# conda_include = []
# conda_lib = []
# if conda_path:
#     conda_base = conda_path.replace("bin/conda", "")
#     conda_include.append(os.path.join(conda_base, "include"))
#     conda_lib.append(os.path.join(conda_base, "lib"))
# lib_dirs.extend(conda_lib)
# include_dirs.extend(conda_include)

extensions = [
    Extension(
        "bamread.src.bamread",
        ["bamread/src/bamread.pyx"],
        # language="c++",
        include_dirs=include_dirs,
        library_dirs=lib_dirs,
        extra_compile_args=compile_options)
]
# libraries=["z"])]

setup(
    name="bamread",
    packages=find_packages(),
    ext_modules=cythonize(extensions, annotate=True, language_level='3'),
    scripts=["bin/bamread"],
    package_data={'': ['*.pyx', '*.pxd', '*.h', '*.c']},
    version=__version__,
    description="Read bam files quickly into dataframes in Python",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/bamread",
    license=["MIT"],
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta", "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    include_dirs=["."],
    long_description=open("README.md").read())
