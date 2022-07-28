## Notes Uniprot API
source: [EMBL-EBI Uniprot API video guide (2021)](https://embl-ebi.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=efde16e9-c2a9-427a-ada2-adce007a9df1)
[slides](https://drive.google.com/file/d/1qZXLl19apibjszCXroC1Jg2BGjcychqX/view)

Data sources:  
FTP, big one-off download.  
**API, medium-size download.**  
Website download, manual small one-off downloads.  

New backend for new website 2021
Anything on the website can be downloaded as an end-point

**Uniprot has various endpoints**, the main data being:
- UniProtKB ()
- UniRef (clusters)
- niParc (sequence archives)

**Supporting data:**
- Taxonomy
- Keywords
- Citations
- Diseases
- Cross-references
- Subcellular locations

**Tools**  
- BLAST
- Peptide search

## Documentation for the new Uniprot API
[OpenAPI - Uniprot](https://www.uniprot.org/help/api)

uniprotkb:  
**GET** /uniprotkb/search

Searh query format: Use advanced search on website
