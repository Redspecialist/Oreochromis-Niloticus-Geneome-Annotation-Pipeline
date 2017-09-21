The following pipeline was developed as a means to correlate all LOCIDs collected from runs of WGCNA using genes from the Oreochromis niloticus genome with their corresponding transcript accession ids, Linkage Groups, and GO terms. Once a reference file is made mapping relevant GOTerms to each gene, filtration is done using a library file that contains different weighting for each GOTerm.
There is added functionality for comparing two networks and creating intersections between the modules of two different WGCNA runs. Intersections with important genes, like the opsins and relevant transcription factors, are flagged by this tool.

The input to this pipeline is obtained from the output of WGCNA and should be input aas a tsv with the following layout. Quotes around each piece of data are handeled and not required:
1 LINE HEADER
"<Gene Accsession ID>"	   "<module number>"	  "<module color>"

Usage of each pipline tool:

go_module_restriction.sh <WGCNA output 1> <WGCNA output 2> <output prefix>

Creates two reference files and GO analysis files for WGCNA outputs. The two reference files in both their raw and Go-term filtered forms are then compared to determine the intersection between two modules.
The resulting matrix is <WGCNA output 1> along the top and <WGCNA output 2> along the side.
References that already exist are not recreated.

single_file_module_interpretation.sh <WGCNA output>

Creates a GO-reference file for one WGCNA output.
