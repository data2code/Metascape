#!/usr/bin/env python
import shlex, subprocess
import numpy as np
import pandas as pd
import math
import os
import shutil
import re
import sys
import util
import urllib
from os import sys, path
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../'))    

from IPython.core.debugger import Tracer
import cPickle

class MyUtils(object):


    @staticmethod
    def dump_object(obj, s_name='untitled', s_cache_dir="."):
        """use to dump any object
        s_name: str, filename, defaults to 'untitled', filename will be appended with suffix '.pickle'
        s_cache_dir: str, filepath, defaults to current folder '.'
        no return"""
        if not s_name.endswith('.pickle'):
            s_name+='.pickle'
        import os
        f=open(os.path.join(s_cache_dir, s_name), 'w')
        cPickle.dump(obj, f)
        f.close()

    @staticmethod
    def load(s_file):
        """load a python object from a pickle file (resulted from a previous dump.
        s_file: str, full pickle file name
        return object
        Example: HCIObject.load('cache/mydump.pickle')"""
        if not s_file.endswith('.pickle'):
            s_file+='.pickle'
        f=open(s_file)
        x=cPickle.load(f)
        f.close()
        return x

