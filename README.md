# This is a Gene Set Enrichment Analysis performed on a ranked list of genes which are differentially expressed in acute lymphocytic leukemia (ALL) cells compared to acute myeloid leukemia (AML) cells. 

The dataset contains 48 samples of patient data, 24 of which are ALL cells and 24 are AML cells. 

The gene expression levels measured in the 48 samples are ranked by their p-values multiplied by the level of up-or downregulation (log 2-fold change). In other words, when a particular gene is significantly downregulated in ALL cells, it is expressed more often in AML cells than in ALL cells.

Download the files below:
* *Leukemia_collapsed_symbols.gct.txt*, file containing the genes
* *Leukemia.cls.txt*, file containing the phenotypes  
* *c2.all.v7.4.symbols.gmt*, gene set file

If you are working in Google Colab, make sure you have uploaded the files to the files tab. If you are working locally, you can run line 65 if the files are located within the same folder as this notebook. Otherwise, add the path that leads to the folder in which you have stored the files (e.g. "C:/Downloads/") in `path_to_files = "C:/Downloads/"`.
