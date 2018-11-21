#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re
from six.moves import zip

class FetchList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'PubmedArticle', xml)

    def attributes(self):
        return ['pubmed_id', 'journal', 'title', 'abstract', 'author']

    def get_uid(self):
        return 'pubmed_id'

    def to_list(self, S_attr=None):
        S_attr=['title', 'journal'] if S_attr is None else S_attr
        if 'pubmed_id' not in S_attr:
            S_attr.append('pubmed_id')
        return el.EntityList.to_list(self, S_attr)

    def pubmed_id(self):
        return self._get_value('./MedlineCitation/PMID')

    def journal(self):
        S_title=self._get_value('./MedlineCitation/Article/Journal/Title')
        S_year=self._get_value('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
        S_month=self._get_value('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Month')
        S_day=self._get_value('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Day')
        S_volume=self._get_value('./MedlineCitation/Article/Journal/JournalIssue/Volume')
        S_issue=self._get_value('./MedlineCitation/Article/Journal/JournalIssue/Issue')
        S_page=self._get_value('./MedlineCitation/Article/Pagination/MedlinePgn')
        return {'title':S_title, 'year':S_year, 'month':S_month, 'day':S_day, 'volume':S_volume, 'issue':S_issue, 'page':S_page}

    def title(self):
        return self._get_value('./MedlineCitation/Article/ArticleTitle')

    def abstract(self):
        return self._get_value('./MedlineCitation/Article/Abstract/AbstractText')

    def author(self):
        out=self._get_nodes('./MedlineCitation/Article/AuthorList/Author')
        S_author=[]
        for x in out:
            S_last=[]
            S_first=[]
            for y in x:
                if y.find('./LastName') is not None:
                    S_last.append(y.find('./LastName').text)
                else:
                    S_last.append('')
                if y.find('./Initials') is not None:
                    S_first.append(y.find('./Initials').text)
                else:
                    S_first.append('')
            S_author.append([x+" "+y for x,y in zip(S_last, S_first)])
        return S_author

class SummaryList(el.EntityList):

    def __init__(self,xml, version=2):
        self.version=version
        if self.version>1:
            el.EntityList.__init__(self, 'DocumentSummarySet/DocumentSummary', xml)
        else:
            el.EntityList.__init__(self, 'DocSum', xml)

    def attributes(self):
        return ['pubmed_id', 'journal', 'title', 'author']

    def get_uid(self):
        return 'pubmed_id'


    def to_list(self, S_attr=None):
        S_attr=['title', 'journal'] if S_attr is None else S_attr
        if 'pubmed_id' not in S_attr:
            S_attr.append('pubmed_id')
        return el.EntityList.to_list(self, S_attr)

    def pubmed_id(self):
        if self.version>1:
            return self._get_attrib('.', 'uid')
        else:
            return self._get_value('./Id')

    def journal(self):
        if self.version>1:
            S_title=self._get_value("./FullJournalName")
            S_date=self._get_value("./PubDate")
        else:
            S_title=self._get_value("./Item[@Name='FullJournalName']")
            S_date=self._get_value("./Item[@Name='PubDate']")
        S_year=[None]*len(S_date)
        S_month=[None]*len(S_date)
        S_day=[None]*len(S_date)
        for i,x in enumerate(S_date):
            if x is None: continue
            y=re.match(r'(?P<year>\d{4})( (?P<month>\D{3})( (?P<day>\d{1,2}))?)?', x)
            S_year[i]=y.group('year')
            S_month[i]=y.group('month')
            S_day[i]=y.group('day')
        if self.version>1:
            S_volume=self._get_value("./Volume")
            S_issue=self._get_value("./Issue")
            S_page=self._get_value("./Pages")
        else:
            S_volume=self._get_value("./Item[@Name='Volume']")
            S_issue=self._get_value("./Item[@Name='Issue']")
            S_page=self._get_value("./Item[@Name='Pages']")
        return {'title':S_title, 'year':S_year, 'month':S_month, 'day':S_day, 'volume':S_volume, 'issue':S_issue, 'page':S_page}

    def title(self):
        if self.version>1:
            return self._get_value("./Title")
        else:
            return self._get_value("./Item[@Name='Title']")

    def author(self):
        if self.version>1:
            out=self._get_nodes("./Authors/Author/Name")
        else:
            out=self._get_nodes("./Item[@Name='AuthorList']/Item[@Name='Author']")
        S_author=[]
        for x in out:
            S_author.append([y.text for y in x])
        return S_author

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.efetch({'db':'pubmed', 'id':'21524757,25100532'})
    x=FetchList(out)
    print(">>>Fetch")
    print(x.to_list(x.attributes()))

    out=eu.esummary({'db':'pubmed', 'id':'21524757,25100532'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))
    exit()
