#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'DocSum', xml)

    def get_uid(self):
        return 'o_id';
    def attributes(self):
        return ['o_id', 'title', 'alt_titles', 'locus']

    def to_list(self, S_attr=None):
        S_attr=['title', 'alt_titles', 'locus'] if S_attr is None else S_attr
        if 'o_id' not in S_attr:
            S_attr.append('o_id')
        return el.EntityList.to_list(self, S_attr)

    def o_id(self):
        return self._get_value('./Id')

    def title(self):
        return self._get_value("./Item[@Name='Title']")

    def alt_titles(self):
        return self._get_values("./Item[@Name='AltTitles']")

    def locus(self):
        return self._get_value("./Item[@Name='Locus']")


FetchList=SummaryList
if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.esummary({'db':'omim', 'id':'219700'})
    print(out)
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))

