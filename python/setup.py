#!/usr/bin/env python

from setuptools import setup

## get long description from README
#with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name='vacuumms',
    version='1.2.1.post2',  
    description='A module for working with vacuumms',
    long_description='A module for working with vacuumms',
    url='http://www.vacuumms.org',
    author='Frank T. Willmore',
    author_email='frankwillmore@gmail.com',
    license='MIT License',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering',
#        'Operating System :: MacOS :: MacOS X',
#        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux'],
    install_requires=['numpy', 'scipy'],
#    extras_require={
#        "docs": ["sphinx", "sphinx_rtd_theme", "numpydoc"]
#    },
#    tests_require=["pytest"],
    packages=[
         'vacuumms',
#        'poo',
#        'poo.poo',
         ],
    include_package_data=True,
)
