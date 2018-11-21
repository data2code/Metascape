#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
import util
from core import *
#from IPython.core.debugger import Tracer
from gputil import GPUtils


class UrlDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.url= xe.attrib['source']
        if not 'dest' in xe.attrib:
            tmp = urlparse.urlsplit(self.url)
            filename, file_ext = splitext(basename(tmp.path))
            self.dest=filename+file_ext
        else:
            self.dest=xe.attrib['dest']

    def get_fn_dest(self):
        return os.path.join(SyncDB.DOWNLOAD_DIR(),self.dest);

    def populate_more(self,root):
        XmlClass.populate_more(self,root)
        self.inputs=['ds:%s' % self.ds]

    def do_update(self):
        #dirty =  self.is_dirty(self.fn_dest)|SyncDB.IS_REBUILD
        dirty =  self.is_dirty(self.fn_dest)
        if dirty:
            urllib.urlretrieve(self.url, self.fn_dest)
        #if SyncDB.IS_REBUILD:
         #   return True
        return dirty

    def get_id(self):
        return self.tag +":" + self.dest
        
    def check_inputs(self):
        return GPUtils.check_url(self.url)
