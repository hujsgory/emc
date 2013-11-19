from distutils.core import setup, Extension
import os

c_ext=Extension("_smn",
                include_dirs = ["c:\\WinPython\\python-2.7.5.amd64\\Lib\\site-packages\\numpy\\core\\include"],
                library_dirs = ["C:\\WinPython\\python-2.7.5.amd64\\Lib\\site-packages\\numpy\\core\\lib"],
                sources      = ["emc/_smn.cpp"])

setup(name         =  'emc'                ,
      version      =  '0.0.2'              ,
      packages     = ['emc']               ,
      package_dir  = {'emc':'emc'}    ,
      package_data = {'emc':['test/*.txt','test/*.py','examples/*']},
      ext_modules  = [c_ext],
      author       =  'hujsgory'           ,
      download_url =  'https://github.com/hujsgory/emc'
      )