#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension
from pkg_resources import resource_filename
from subprocess import call as subpcall
from setuptools import find_packages

if sys.version < "2.6.0" or sys.version > "2.8.0":
    print "Please use a Python with higher version than 2.6.0"
    sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
    sys.exit(1)
    exit(1)
    
def run_cmd(command):
    subpcall (command, shell = True)
def compilemis():
    curdir = os.getcwd()
    os.chdir('BETA/misp')
    run_cmd('make')
    run_cmd('chmod 755 *')
    os.chdir(curdir)
    
def main():

    #compilemis()
    if not float(sys.version[:3])>=2.5:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.5! python 2.6.1 or newer is recommended!\n")
        sys.exit(1)
    setup(name="BETA-Package",
          version="1.0.7",
          description="BETA -- Binding and Expression Targets Analysis ",
          author='Su Wang',
          author_email='wangsu0623@gmail.com',
          package_dir={'BETA' : 'BETA'},
          install_requires=['argparse','numpy'],
          packages=['BETA'],
          scripts=['bin/BETA','BETA/misp/misp'],
          package_data={'BETA':['references/*','templates/*']},

          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Database',
            ],
          )
    
    
if __name__ == '__main__':
    compilemis()
    main()
    
