#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
from core import SyncDB,XmlClass
#from .. import SyncDB
class ColConst(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe=xe)
        self.col_name = xe.attrib['colName']
        self.value = xe.attrib['value']
        self.required =True

    def do_work(self,row=None):
        return self.value
