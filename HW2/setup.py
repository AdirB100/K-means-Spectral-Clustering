from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    version="0.1.0",
    author="Adir Barda 318574894, Offir Dassa 315326652",
    author_email="adirbarda@mail.tau.ac.il , offirdassa@mail.tau.ac.il",
    packages=find_packages(),
    license='GPL-2',
    ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])])

