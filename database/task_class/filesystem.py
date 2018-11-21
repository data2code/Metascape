import numpy as np
import shutil
import urllib
import urlparse
import os.path as path
import os
from core import *
from IPython.core.debugger import Tracer

class FileSystem(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.source = xe.attrib['source']
        if path.dirname(self.source) == '':
            self.source = path.join(SyncDB.DOWNLOAD_DIR(), self.source);
        if not 'dest' in xe.attrib:
            self.dest = path.join(SyncDB.DOWNLOAD_DIR(), path.basename(self.source));
        else:
            self.dest= path.join(SyncDB.DOWNLOAD_DIR(), xe.attrib['dest']);
        
    def get_fn_dest(self):
        return self.dest;

    def populate_more(self,root):
        XmlClass.populate_more(self,root)
        self.inputs=['ds:%s' % self.ds]

    def do_update(self):
        #dirty =  self.is_dirty(self.fn_dest)|SyncDB.IS_REBUILD
        if not path.exists(self.source):
            print 'FileSystem: source file "',  self.source, '" does not exist'
            sys.exit()
        dirty =  self.is_dirty(self.fn_dest)
        if dirty and self.source != self.dest:
            cmd = 'cp ' + str(self.source) + ' ' + str(self.dest);
            util.unix(cmd);
        #if SyncDB.IS_REBUILD:
         #   urn True
        return dirty

    def get_id(self):
        return self.tag +":" + self.dest
    
    def check_inputs(self):
        if not path.exists(self.source):
            print 'FileSystem: source file "',  self.source, '" does not exist'
            return False
        return True
            
