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


class IdMapping(UploadCsvConvert):
    def __init__(self, xe):
        #add default col instance to children
        xe.attrib['newCols'] = 'gid,id_type_id,source_id,id_status,ds,tax_id'
        UploadCsvConvert.__init__(self,xe=xe,dest='gid2source_id')
        self.type_col = 'id_type_id'
        if 'destFile' in xe.attrib:    
            self.use_dest =xe.attrib['destFile'];

    def get_type_col_value_sql(self):
        s = 'SELECT id_type_id FROM %s.id_type WHERE id_type_name = ?' % SyncDB.DATABASE
        return s



