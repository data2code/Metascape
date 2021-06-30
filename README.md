# Metascape
Metascape Code Repository.

For the latest version of the code, please see: https://github.com/data2code/msbio/

# old instructions

The example script below relies on a backend MySQL database to run.

```
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
```

Output

```
MetaEngine> Start timer ...
output:8661|MetaEngine> Passed: 0.0 secs, Membership
output:8661|MetaEngine> Passed: 13.8 secs, genelists_overlap_analysis
Enrichment analysis of individual lists ...
Loading Entrez gene table for 9606...
Loading mapping table for 9606...
GeneLists::go_analysis> Start timer ...
output:8661|MetaEngine> Passed: 16.2 secs, GO membership: t_go_mem 1635 x 11
output:8661|MetaEngine> Passed: 0.0 secs, GO Enrichment
Interactome analysis of individual lists ...
loading PPI database from database for tax_id: 9606 ...
output:8661|MetaEngine> Passed: 14.1 secs, Start PPI Analysis
GeneList::ppi_analysis> Start timer ...
GeneList::ppi_analysis> Start timer ...
GeneList::ppi_analysis> Start timer ...
Find MCODE components ...
Find MCODE components ...
Find MCODE components ...
Found 10 MCODE components
Found 1 MCODE components
output:8661|GeneList::ppi_analysis> Passed: 4.8 secs, Network DONE
Found 7 MCODE components
output:8661|GeneList::ppi_analysis> Passed: 5.2 secs, Network DONE
Found 9 MCODE components
output:8661|GeneList::ppi_analysis> Passed: 7.4 secs, Network DONE
Evidence collecting ...
Meta-analysis, combined enrichment analysis ...
MetaEngine::GPEC_go_analysis> Start timer ...
GO_Cluster> Start timer ...
output:8661|GO_Cluster> Passed: 2.8 secs, Done membership...
GO_Cluster::cluster> Start timer ...
Matrix size: 709 genes x 1555 GOs
output:8661|GO_Cluster::cluster> Passed: 2.4 secs, Kappa done ...
GOList::membership> Start timer ...
output:8661|GOList::membership> Passed: 0.1 secs, Score calculation
output:8661|GOList::membership> Passed: 0.0 secs, Get membership matrix
output:8661|GOList::membership> Passed: 1.7 secs, Cluster genes and columns for each cluster.
output:8661|MetaEngine> Passed: 31.5 secs, Circos GO Overlap
Meta_analysis, Combined interactome analysis ...
Use GPEC Network> GPEC size=426, MERGE size=845
Find MCODE components ...
Found 13 MCODE components
output:8661|MetaEngine::GPEC_ppi_analysis> Passed: 18.8 secs, ***Start GO Analysis***
> Start timer ...
output:8661|MetaEngine::GPEC_ppi_analysis> Passed: 4.0 secs, Finish GO analysis
output:8661|MetaEngine::GPEC_ppi_analysis> Passed: 0.6 secs, *****GO Analysis within PPI Analysis
Check the output folder: output
```
