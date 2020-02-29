In this folder, we have executable workflows (.sh or .R) and R notebooks. You should execute them in the order they are numbered.

0: This loads files and filters OTU tables, metadata files, fasta files, etc. 
1: This step clusters OTUs against the woodhams anti-bd database at 100% identity using vsearch
2: This step adds microbial community parameters to the metadata and adjusts parameter values so they are compliant with bayesian analysis
3: This is an exploratory data analysis step that confirms the assumptions and distribution of data for each microbial community parameter
4: This runs the bayesian analysis that creates the posterior distributions for each microbial community metric
5: This runs the statistics for all analyses
5: This creates the figures for manuscript

26Nov2019
- Currently re-running SILVA picking of OTUs. After this, should re-cluster woodhams database. 
- In the meantime, workng on #2 workflow.
