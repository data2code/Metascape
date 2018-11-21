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
from pprint import pprint
import StringIO
import db
from gp import *
from core import *
from IPython.core.debugger import Tracer

class Interaction(UploadCsvConvert):
    def __init__(self, xe):
        xe.attrib['newCols'] = 'ds,gid_A,gid_B,interaction_category,interaction_categy,interaction_type,interaction_type_id,score,support,tax_id_A,tax_id_B'
        UploadCsvConvert.__init__(self,xe=xe,dest='interaction')
        self.type_col = 'interaction_type_id'

    def get_type_col_value_sql(self):
        return 'SELECT interaction_type_id FROM %s.interaction_type WHERE interaction_type_name = ?' % SyncDB.DATABASE



