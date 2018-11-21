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


class GeneGoTerms(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe)
        type_names = [ x.strip() for x in xe.attrib['typeNames'].split(',')]
        for type_name in type_names:
            term = """<Term source="genego_terms.csv" 
                            typeName='%s' 
                            filter_str='type_name in(&quot;%s&quot;)' >
                      </Term>"""            
            term_gene_pair = """<TermGenePair source='genego_term_gene_pair.csv'    
                                    typeName='%s' 
                                    filter_str='type_name in(&quot;%s&quot;)' />"""

            term_xml = term % (type_name,type_name)
            term_root = ET.fromstring(term_xml)
            term_c = SyncDB.task_factory(term_root)
            self.children.append(term_c)
            term_gene_pair_xml = term_gene_pair % (type_name,type_name)
            term_gene_pair_root = ET.fromstring(term_gene_pair_xml)
            term_gene_pair_c = SyncDB.task_factory(term_gene_pair_root)
            self.children.extend(term_gene_pair_c.get_xmlclass_instances())

    def get_xmlclass_instances(self):
        r = list(self.children)
        return r
