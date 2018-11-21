#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re

# mesh does not return XML format for efetch

class SummaryList(el.EntityList):

    def __init__(self,xml, version=2):
        self.version=version
        if self.version>1:
            el.EntityList.__init__(self, 'DocumentSummarySet/DocumentSummary', xml)
        else:
            el.EntityList.__init__(self, 'DocSum', xml)

    def get_uid(self):
        return 'mesh_id'

    def attributes(self):
        return ['mesh_id', 'mesh_ui', 'year', 'terms', 'description', 'parent', 'children','tree_num']

    def to_list(self, S_attr=None):
        S_attr=['terms','description','tree_num'] if S_attr is None else S_attr
        if 'mesh_id' not in S_attr:
            S_attr.append('mesh_id')
        if 'mesh_ui' not in S_attr:
            S_attr.append('mesh_ui')
        return el.EntityList.to_list(self, S_attr)

    def mesh_id(self):
        if self.version>1:
            return self._get_attrib('.', 'uid')
        else:
            return self._get_value('./Id')

    def mesh_ui(self):
        if self.version>1:
            return self._get_value('./DS_MeSHUI')
        else:
            return self._get_value('./Id')

    def year(self):
        if self.version>1:
            return self._get_value("./DS_YearIntroduced")
        else:
            return self._get_value("./Item[@Name='DS_YearIntroduced']")

    def terms(self):
        if self.version>1:
            return self._get_value("./DS_MeSHTerms/string")
        else:
            return self._get_values("./Item[@Name='DS_MeshTerms']/Item[@Name='string']")

    def description(self):
        if self.version>1:
            return self._get_value("./DS_ScopeNote")
        else:
            return self._get_value("./Item[@Name='DS_ScopeNote']")

    def parent(self):
        if self.version>1:
            return self._get_value("./DS_IdxLinks/LinksType/Parent")
        else:
            return self._get_value("./Item[@Name='DS_IdxLinks']/Item[@Name='LinksType']/Item[@Name='Parent']")

    def children(self):
        if self.version>1:
            return self._get_values("./DS_IdxLinks/LinksType/Children/int")
        else:
            return self._get_values("./Item[@Name='DS_IdxLinks']/Item[@Name='LinksType']/Item[@Name='Children']/Item[@Name='int']")

    def tree_num(self):
        if self.version>1:
            return self._get_value("./DS_IdxLinks/LinksType/TreeNum")
        else:
            return self._get_value("./Item[@Name='DS_IdxLinks']/Item[@Name='LinksType']/Item[@Name='TreeNum']")





if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()

    out=eu.esummary({'db':'mesh', 'id':'68017209,68017220'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))

