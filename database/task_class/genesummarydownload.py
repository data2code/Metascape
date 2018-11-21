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
from core import *
import util
import csv
import urllib2
import json
import sys
from IPython.core.debugger import Tracer
#tew

class GeneSummaryDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_source = os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['source']); 
        self.fn_dest =os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['dest']) 
        self.inputs = [self.fn_source]
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;

    def populate_more(self,root):
        self.outputs = [self.fn_dest]           
        
    def do_update(self):
        dirty = True;
        import os.path, time,datetime
        if os.path.isfile(self.fn_dest): 
            statinfo = os.stat(self.fn_dest)
            if statinfo.st_size > 0:
                file_time = os.path.getmtime(self.fn_dest)
                current_time = datetime.datetime.now()
                time_delta_hour = (current_time -datetime.datetime.fromtimestamp(file_time)).total_seconds()/3600
                print 'time_delta_hour for gene info file', time_delta_hour
                if time_delta_hour <= 100*24:
                    dirty  = False

        if (dirty):            
            from eutils_summary import EUtilsSummary
            EUtilsSummary.get_gene_info(self.fn_source, self.fn_dest, self.taxidList)

