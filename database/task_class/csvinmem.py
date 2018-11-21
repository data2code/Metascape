#!/usr/bin/env python
#from .core import *
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
import util
from pprint import pprint
import StringIO
import db
from gp import *
from core import *


class CsvInMem(CsvConvert):
    def __init__(self, xe):
        CsvConvert.__init__(self,xe=xe)
        self.dest = self.__class__.__name__+"_"+xe.attrib['source'].replace('.gz','')

    
    def populate_more(self,root):
        CsvConvert.populate_more(self,root)
        self.outputs=[self.dest,self.fn_dest]

    def get_fn_dest(self):
        return os.path.join(SyncDB.DOWNLOAD_DIR(),self.dest)



