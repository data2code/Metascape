#!/usr/bin/env python
import pandas as pd
import ms.biolists as bl
import ms.meta_engine as me
import os
import go

# create three influenza hit lists
t=pd.read_csv('influenza_hits.csv')
lists=bl.GeneLists([ bl.GeneList(t.ix[i, 'Name'], t.ix[i,'Hits'].split(' ')) for i in t.index ])

# define output folder
s_output="output"
if not os.path.exists(s_output):
    os.mkdir(s_output)

# create an analysis engine
eng=me.MetaEngine(s_output)

# specify the GeneLists object
eng.puts({'lists':lists})

# analyze gene membership across all hit lists
eng.genelists_overlap_analysis()

print("Enrichment analysis of individual lists ...")
eng.go_analysis()
print("Interactome analysis of individual lists ...")
eng.ppi_analysis()
print("Evidence collecting ...")
eng.GPEC_evidence()

print("Meta-analysis, combined enrichment analysis ...")
eng.GPEC_go_analysis()
print("Meta_analysis, Combined interactome analysis ...")
eng.GPEC_ppi_analysis()

print("Check the output folder: %s" % s_output)
