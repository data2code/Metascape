#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
import re
from col import Col

from core import *

class  ColNanStr(Col):
    def __init__(self, xe):
        Col.__init__(self,xe=xe)

    def do_work(self,row=None):
        r=Col.do_work(self,row)
        r= str(r).replace('.0','')
        return r
