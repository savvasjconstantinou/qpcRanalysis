script to analyze qPCR data using deltadeltaCq method

Data needs to be read in as csv file. To create csv file: On a new spreadsheet reorganize the data to eliminate all spaces, unused wells and preliminary notes produced by thermocycler. Columns will also need to be added, which at minimum should include: "Gene", "GeneID", "Sample", "SampleID", and "Cq". Gene column should indicate name of gene examined (ex: scn4aa/Bactin etc.). GeneID should be denoted "GOI" for our gene of interest (scn4aa in this case) and "Ref1" for Bactin and "Ref2" for rsp11 (Ref3 etc...if additional). Sample column should indicate which specimen is being utilized (ex: MOC/MOL etc.). SampleID should be denoted "control", "exp1" or "exp2"- example: MOC is the control, MOL and MOH would be experimental.

Below is an example of what the final edited .csv file could look like:

Well	Cq	Gene	GeneID	Sample 	SampleID
A04	31.64896	scn4aa	GOI	other	exp1
A06	32.31516	scn4aa	GOI	other	exp1
B04	24.33583	scn4aa	GOI	MOC	control
B05	25.06438	scn4aa	GOI	MOC	control
B06	22.79898	scn4aa	GOI	MOC	control
C04	26.22798	bactin	Ref1	other	exp1
C05	25.29998	bactin	Ref1	other	exp1
C06	25.47763	bactin	Ref1	other	exp1
D04	22.14269	bactin	Ref1	MOC	control
D05	21.6308	bactin Ref1	MOC	control
D06	20.82959 bactin	Ref1	MOC	control

The graph will need to be changed to have correct titles and call particular variables of interest on a case by case instance.
