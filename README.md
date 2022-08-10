# presconvis
**Work in progress**  
Quick and simple PRotein RESidue CONservation and VISualization tool.

<a href="https://github.com/DavidL-H/presconvis/pulse" alt="Activity"><img src="https://img.shields.io/github/commit-activity/m/DavidL-H/presconvis" /></a>
![Python 3](https://img.shields.io/badge/Language-Python_3-red.svg)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)

This tool takes a single protein sequence and user input, grabbing n proteins from uniprot (UnirefXX), performing MSA, and visualizing both a tree for total conservation, but also a 2d sequence plot, showing indivudual residue conservation.

The front end will be build in with the newly launched Shiny for Python.
Persistent data storage wil MySQL

This tools accesses protein sequences and their meta data from the **Uniprot REST API**:  
https://rest.uniprot.org/  
https://www.uniprot.org/help/api  

Multiple sequence alignments are run remotely with Clustal Omega via the **EMBL-EBI REST API**:  
https://www.ebi.ac.uk/Tools/services/rest/clustalo  
https://www.ebi.ac.uk/Tools/common/tools/help/  
