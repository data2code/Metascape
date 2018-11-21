#!/usr/bin/env python
#from .core import *
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from core import *

class ColUnique(Col):
    def __init__(self, xe):
        Col.__init__(self,xe=xe)

    def do_work(self,row=None):
        s = row[self.source_col]
        """ :type s:str"""
        vs = '; '.join({ x.strip() for x in   s.split(r';')})
        return vs
