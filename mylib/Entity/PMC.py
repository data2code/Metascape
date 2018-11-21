#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re
from six.moves import zip

class FetchList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'article', xml)

    def attributes(self):
        return ['pmc_id', 'pubmed_id', 'journal', 'title', 'abstract', 'author']

    def to_list(self, S_attr=None):
        S_attr=['title', 'journal'] if S_attr is None else S_attr
        if 'pmc_id' not in S_attr:
            S_attr.append('pmc_id')
        return el.EntityList.to_list(self, S_attr)

    def pmc_id(self):
        return self._get_value("./front/article-meta/article-id[@pub-id-type='pmc']")

    def pubmed_id(self):
        return self._get_value("./front/article-meta/article-id[@pub-id-type='pmid']")

    def journal(self):
        S_title=self._get_value('./front/journal-meta/journal-title-group/journal-title')
        S_year=self._get_value('./front/article-meta/pub-date/year')
        S_month=self._get_value('./front/article-meta/pub-date/month')
        S_day=self._get_value('./front/article-meta/pub-date/day')
        S_volume=self._get_value('./front/article-meta/volume')
        S_issue=self._get_value('./front/article-meta/issue')
        S_page=self._get_value('./front/article-meta/fpage')
        S_lpage=self._get_value('./front/article-meta/lpage')
        for i,x in enumerate(S_page):
            if x is not None and S_lpage[i] is not None:
                S_page[i]+='-'+S_lpage[i]
        return {'title':S_title, 'year':S_year, 'month':S_month, 'day':S_day, 'volume':S_volume, 'issue':S_issue, 'page':S_page}

    def title(self):
        return self._get_value('./front/article-meta/title-group/article-title')

    def abstract(self):
        return self._get_value('./front/article-meta/abstract/p')

    def author(self):
        out=self._get_nodes("./front/article-meta/contrib-group/contrib[@contrib-type='author']/name")
        S_author=[]
        for x in out:
            S_last=[]
            S_first=[]
            for y in x:
                if y.find('./surname') is not None:
                    S_last.append(y.find('./surname').text)
                else:
                    S_last.append('')
                if y.find('./given-names') is not None:
                    S_first.append(y.find('./given-names').text)
                else:
                    S_first.append('')
            S_author.append([x+" "+y for x,y in zip(S_last, S_first)])
        return S_author

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'DocSum', xml)

    def attributes(self):
        return ['pmc_id', 'pubmed_id', 'journal', 'title', 'author']

    def to_list(self, S_attr=None):
        S_attr=['title', 'journal'] if S_attr is None else S_attr
        if 'pmc_id' not in S_attr:
            S_attr.append('pmc_id')
        return el.EntityList.to_list(self, S_attr)

    def pmc_id(self):
        return self._get_value('./Id')

    def pubmed_id(self):
        return self._get_value("./Item[@Name='ArticleIds']/Item[@Name='pmid']")

    def journal(self):
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
        S_volume=self._get_value("./Item[@Name='Volume']")
        S_issue=self._get_value("./Item[@Name='Issue']")
        S_page=self._get_value("./Item[@Name='Pages']")
        return {'title':S_title, 'year':S_year, 'month':S_month, 'day':S_day, 'volume':S_volume, 'issue':S_issue, 'page':S_page}

    def title(self):
        return self._get_value("./Item[@Name='Title']")

    def author(self):
        out=self._get_nodes("./Item[@Name='AuthorList']/Item[@Name='Author']")
        S_author=[]
        for x in out:
            S_author.append([y.text for y in x])
        return S_author

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.efetch({'db':'pmc', 'id':'4160113,4160114'})
    x=FetchList(out)
    print(">>>Fetch")
    print(x.to_list(x.attributes()))

    out=eu.esummary({'db':'pmc', 'id':'4160113,4160114'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))
    exit()

