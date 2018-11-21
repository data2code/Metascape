#!\usr\bin\env python

import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
from core import *
import util
import csv
from gputil import GPUtils
from lxml import etree

from IPython.core.debugger import Tracer

class BioconductorAffymetrix(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.tag = "BioconductorAffymetrix"
        self.fn_dest=os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['dest'])
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest]
   
    def do_update(self):
        cmd = "time " + os.path.dirname(os.path.realpath(__file__)) + "/Bioconductor/affy_array.R " + self.fn_dest;
        print cmd;
        util.unix(cmd);
                       
    def check_inputs(self):
        return True;