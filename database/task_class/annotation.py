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

class Annotation(UploadCsvConvert):
    def __init__(self, xe):
        xe.attrib['newCols'] = 'gid,annotation_type_id,content,annotation_field1,ds,tax_id'
        UploadCsvConvert.__init__(self,xe=xe,dest='annotation')
        self.type_col = 'annotation_type_id'

    def get_type_col_value_sql(self):
        return 'SELECT annotation_type_id FROM %s.annotation_type WHERE annotation_type_name = ?' % SyncDB.DATABASE



