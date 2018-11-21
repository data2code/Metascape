#!/usr/bin/env python
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
import util
from pprint import pprint
import StringIO
import db
from core import *


class Main(XmlClass):
    def __init__(self, xe):
        #add default col instance to children
        XmlClass.__init__(self,xe=xe)

