#!/usr/bin/env python
from os import sys, path
p1 = path.join(path.dirname(path.abspath(__file__)),'../mylib')
print p1
sys.path.insert(0, p1)
import pandas as pd
import util
import db
import urllib


con=db.get_con('GENEGO')
t = db.from_sql('select 1 from dual')
print t
