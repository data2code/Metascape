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

class Homologene(UploadCsvConvert):
    def __init__(self, xe):
        UploadCsvConvert.__init__(self,xe=xe,dest='homologene')



