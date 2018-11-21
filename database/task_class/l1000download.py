#!/usr/bin/env python
import numpy as np
import shutil
import urllib
import urlparse
import os
from core import *
import util
from pprint import pprint
import pandas as pd

class L1000Download(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)

    def populate_more(self,root):
        self.outputs.extend([os.path.join(SyncDB.DOWNLOAD_DIR(), 'L1000_term.csv'),
                             os.path.join(SyncDB.DOWNLOAD_DIR(), 'L1000_term_pair.csv'),
                             ])

    def do_update(self):
        pass
