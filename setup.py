from distutils.core import setup, Extension
import os
from config import *

c_ext=Extension("_smn",
                include_dirs = [NUMPY_ROOT + "core\\include"],
                library_dirs = [NUMPY_ROOT + "core\\lib"],
                sources      = ["emc/_smn.cpp", "emc/_iterative.cpp"])

setup(name         =  'emc'                ,
      version      =  '0.0.2'              ,
      packages     = ['emc']               ,
      package_dir  = {'emc':'emc'}    ,
      package_data = {'emc':['test/*.txt','test/*.py','examples/*']},
      ext_modules  = [c_ext],
      author       =  'hujsgory'           ,
      download_url =  'https://github.com/hujsgory/emc'
      )