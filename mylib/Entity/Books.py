#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re

# books does not return useful XML format for efetch
# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=books&id=1690248

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'DocSum', xml)

    def attributes(self):
        return ['books_id', 'date', 'title', 'type', 'reference_id']

    def to_list(self, S_attr=None):
        S_attr=['title'] if S_attr is None else S_attr
        if 'books_id' not in S_attr:
            S_attr.append('books_id')
        return el.EntityList.to_list(self, S_attr)

    def books_id(self):
        return self._get_value('./Id')

    def date(self):
        S_date=self._get_value("./Item[@Name='PubDate']")
        S_year=[None]*len(S_date)
        S_month=[None]*len(S_date)
        S_day=[None]*len(S_date)
        for i,x in enumerate(S_date):
            if x is None: continue
            y=re.match(r'(?P<year>\d{4})(/0?(?P<month>\d+)(/0?(?P<day>\d{1,2}))?)?', x)
            S_year[i]=y.group('year')
            S_month[i]=y.group('month')
            S_day[i]=y.group('day')
        return {'year':S_year, 'month':S_month, 'day':S_day}

    def title(self):
        return self._get_value("./Item[@Name='Title']")

    def type(self):
        return self._get_value("./Item[@Name='RType']")

    def reference_id(self):
        return self._get_value("./Item[@Name='RID']")

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.esummary({'db':'books', 'id':'1895788,1502792'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))

