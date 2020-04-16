# GLIPH-on-VDJb
###### 06-05-2020
Introduction
----
This repository contains the code and database used to create the gliph results for the research project of S.J. van den Brink .

Database
----
The database has been downloaded from: https://vdjdb.cdr3.net/search on 23-03-2020 using the following filters:
- CDR3 :
  - Species :   
    - Human
  - Genechain : 
    - TRA
    - TRB
    - paired only
    - append pairs
- MHC :
  - Class : 
	  - MHCI
    - MHCII
- Meta :
  - Assay type :  
	  - Multimer sorting
    - Culture-based
    - Other
  - Sequencing :  
	  - Sanger
    - High-throughput
    - Single-cell

Gliph
----
The gliph (Grouping of Lymphocyte Interactions by Paratope Hotspots) algortim has been obtained from: https://github.com/immunoengineer/gliph and the authors article on gliph and its usage can be found here:  doi:10.1038/nature22976

Code
----
The code for the datamanipulation can be found in https://github.com/servan98/GLIPH-on-VDJb/blob/master/Data-manipulation.R \
The code for running gliph can be found on https://github.com/immunoengineer/gliph \
The code for visualisation can be found in https://github.com/servan98/GLIPH-on-VDJb/blob/master/Visualisation_HLAtype.py | https://github.com/servan98/GLIPH-on-VDJb/blob/master/Visualisation_PatientID.py | https://github.com/servan98/GLIPH-on-VDJb/blob/master/Visualisation_Epitope.py | https://github.com/servan98/GLIPH-on-VDJb/blob/master/Visualisation_JgeneSegment.py | https://github.com/servan98/GLIPH-on-VDJb/blob/master/Visualisation_VgeneSegment.py

Packages used
----
In R (Rstudio 1.1.456) : 
- tidyverse
- stringr

In Python (3.7.6) :
- pyvis
- pandas

