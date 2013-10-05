#setup.py
from distutils.core import setup
import os
#data_files=map(lambda y: 'emc/test/'+y,filter(lambda x: '.txt' in x, os.listdir('emc/test')))

setup(name         =  'emc'                ,
      version      =  '0.0.1'              ,
      packages     = ['emc']               ,
      package_dir  = {'emc':'emc'}    ,
      package_data = {'emc':['test/*.txt','test/*.py','examples/*']},
      author       =  'hujsgory'           ,
      download_url =  'https://github.com/hujsgory/emc'
      )