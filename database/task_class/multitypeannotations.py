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
from core import *
#from core import *
from pprint import pprint
import StringIO
import db
from IPython.core.debugger import Tracer

class MultiTypeAnnotations(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe)
        type_names = [ x.strip() for x in xe.attrib['typeNames'].split(',')]
        
        self.tag=xe.attrib['name']
        self.fn_source=xe.attrib['source']
        self.inputs=[self.fn_source]
        for type_name in type_names:
            annotations = """<Annotation source="%s" 
                            typeName='%s' 
                            filter_str='type_name in(&quot;%s&quot;)' >
                      </Annotation>"""            

            annotation_xml = annotations % (self.fn_source, type_name.replace(' ', '_'),type_name) #type_name.replace(' ', '_') is because annotation_type_name cannot have space in it.
            annotation_root = ET.fromstring(annotation_xml)
            annotation_c = SyncDB.task_factory(annotation_root)
            self.children.append(annotation_c)

    def get_xmlclass_instances(self):
        r = list(self.children)
        return r
