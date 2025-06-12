Step-by-step Protocol to Use PepSIFinder
**1. Prepare the Input Data**
Export your peptide-spectrum match (PSM) results from PEAKS or another search engine. Ensure your DataFrame includes at least the following columns: `Peptide`, `RT`, `z`, `1/k0 Range`, `m/z`, `Mass`, `Scan`, and `Accession`.
**2. Load the Data into a DataFrame**
```python
import pandas as pd
df = pd.read_csv("your_psm_data.csv")
```
**3. Run the Candidate Finder**
```python
from PepSIFinder import get_candidates
candidates = get_candidates(df, RT_tol=1.0, ims_diff=0.04)
```
**4. Add Protein Descriptions (Optional)**
```python
from PepSIFinder import add_descriptor
protein_df = pd.read_csv("protein_info.csv")
annotated = add_descriptor(candidates, protein_df)
```
**5. Compare Between Samples (Optional)**
```python
from PepSIFinder import compare_peptide_df
common, only_in_A, only_in_B = compare_peptide_df(df_A, df_B)
```
**7. Find posttranslationally modified Peptides (Optional)**
E.g for finding methylated residues. 
```python
from PepSIFinder import extract_modified_peptides
modified = extract_modified_peptides("your_psm_data.csv", modifications=["+14.02"], residue="D")
```
