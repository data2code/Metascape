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
from core import *


class Term2Term(UploadCsvConvert):
    def __init__(self, xe):
        xe.attrib['newCols'] = 'term_id,parent_term_id, distance,term_category_id' 
        UploadCsvConvert.__init__(self,xe=xe,dest='term2term')
        self.type_col = 'term_category_id'


    def get_type_col_value_sql(self):
        return 'SELECT t.term_category_id FROM %s.term_category t WHERE t.category_name = ?' % SyncDB.DATABASE

