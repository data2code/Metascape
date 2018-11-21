#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
import re
from gp import *
from core import *
import math
from IPython.core.debugger import Tracer

class Col(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe=xe)
        self.col_name = xe.attrib['colName']
        self.source_col = xe.attrib['sourceCol']
        if not 'required' in xe.attrib:
            self.required = False
        else:
            self.required = True
        self.match = None
        if 'match' in xe.attrib:
            self.match = xe.attrib['match']

    def do_work(self,row=None):
        r = str(row.get(self.source_col,None))
        if self.match is not None:
            m = re.search(self.match,r if (r.lower()!='nan')   else '')
            if m and len(m.groups())!=0:               
                r = m.group(1)
            else:
                r = None
        return r
