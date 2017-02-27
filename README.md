# qpcRanalysis
script to analyze qPCR data using deltadeltaCq method

data needs to be read in as csv file.
To create csv file: On a new spreadsheet reorganize the data to eliminate all spaces and add a few columns. At minimum you should include 
columns named, "Gene", "GeneID", "Sample", "SampleID", and "Cq". 
GeneID should be denoted "GOI" for our gene of interest (scn4aa in this case) and "Ref1" for Bactin and "Ref2" for rsp11 (Ref3 etc...if additional)
SampleID should be denoted "control", "exp1" or "exp2"- example: MOC is the control, MOL and MOH would be experimental.

The graph will need to be changed to have correct titles.
