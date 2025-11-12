# NeuroCHARGE_obesity_ctxVol
Scripts to run adiposity-ctx structure associations
# File Instructions

## File 1: `neuro` file

**Notes:** - The column names **and their order** must exactly match those in
`neuro_columns.txt`, which contains names for the FreeSurfer-parcellated regions, left and right hemispheres, and ICV.
- Include the following 212 variables:
1.  `ID` Participant ID (n=1)
2.  Cortical volume for the 68 FreeSurfer-parcellated regions in the left and right hemispheres (n=68)
3.  Cortical thickness for the 68 FreeSurfer-parcellated regions in the left and right hemispheres (n=68)
4.  Cortical surface area for the 68 FreeSurfer-parcellated regions in the left and right hemispheres (n=68)
5.  Total cortical volume for the left and right hemispheres (n=2)
6.  Average cortical thickness in the left and right hemispheres (n=2)
7.  Total cortical surface area for the left and right hemispheres (n=2)
8.  Intracranial volume (ICV) (n=1)

------------------------------------------------------------------------

## File 2: `non.neuro` file

This file must include all non-brain types of variables (i.e., adiposity, genotypes, genetic PCs, any cohort-specific variables).
(See `non_neuro_example.xlsx` for reference.)

| Variable | Description (unit or coding)      | Note                                            |
|-------------|---------------------------------------|---------------------------------------------|
| `ID`     | Participant ID                    |                                                 |
| `FID`    | Family ID (same as `ID` if unrelated)|                                              |
| `age_mri`| Age at MRI scan (years)         |                                                 |
| `age_adiposity`| Age at adiposity measurement (years)|If it was taken at the time of brain MRI, put `age_mri`.|
| `sex`    | Sex (`M/F`)                       |                                                 |
| `icv`    | Intracranial volume (mm³).      |                                                 |
| `current_smoking`| Current smoking staaus (`Yes/No`)|                                          |
| `hypertension`| Has hypertension? (`Yes/No`) |                                                 |
| `T2D`    | Has type 2 diabetes? (`Yes/No`)   |                                                 |
| `BMI`    | Body mass index (kg/m²)           |                                                 |
| `waist`  | Waist circumference (cm)          |                     |
| `height` | Height (cm)                       |                     |
| `WHR`    | Waist Hip Ratio: Waist circumference (cm) / Hip circumference (cm)|                     |
| `Ancestry` | AFR (African), AMR (American), EAS (East Asian), EUR (European), SAS (African)| Contact us if you have non-European ancestries and need to genetically infer ancestries: We will provide scripts. |
| `Ethnicity` | e.g., White, African American, Chinese, ... |                |
| `APOE` | APOE genotypes: coded as `e2e2`, `e2e3`, `e2e4`, `e3e3`, `e3e4`, or `e4e4` |                           |
| `genoPC1` – `genoPC4` or more | The first >=4 leading principal component scores.|At least 4 top PCs will be included. |
| `cohort_specific_1` | e.g., genotype array, MRI site |                     |
| …                   | …                              |                     |
| `cohort_specific_m` | …                              |                     |

------------------------------------------------------------------------

**End of File Instructions**
