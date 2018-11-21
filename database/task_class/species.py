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

#This class is designed to give each tag a customized set of supported species (the global supported species are provieded by DefaultSpecies tag), however b/c we expanded the input tax_id's to include species in taxonomy subtree it not longer works as is right now. Modifications are needed to provide support with species subtree.

class Species(object):
    def __init__(self, xe):
        if 'taxidList' in xe.attrib:
            self.supported_species = filter(None, [x for x in xe.attrib['taxidList'].split(',')])
        else:
            self.supported_species = SyncDB.SUPPORTED_SPECIES
        
        if 'taxidCol' in xe.attrib:
            self.taxid_col = xe.attrib['taxidCol']
        else:
            self.taxid_col = 'tax_id'

	#This function filters out rows that are not of the desired species. However it is no longer used because of the poor peformance. Filtering agains tax_id are done at higher level for better performance.
    def check(self,row=None):
        #Supporting all species.
        if (len(self.supported_species) == 0):
            return True
            
        taxid = str(row.get(self.taxid_col,None))
        if taxid is None:
            return True
        
        for t in taxid.split(','):
            if t in self.supported_species:
                return True
                
        return False
