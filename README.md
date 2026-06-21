# AncestrySkew

Compute ancestry skew for variants in an annotated TSV/TSV.GZ file. The workflow filters variants by PIP, recalculates the GVS max subpopulation from minor allele frequency (MAF), and runs Fisher tests comparing the max subpopulation to the remaining cohort.

The script also reports a second ancestry skew calculation with admixed samples removed from the cohort background. By default, the admixed subpopulation is `oth`, but this can be changed with a comma-separated input such as `oth,amr`.

## Inputs

### WDL workflow

`workflows/ComputeAncestrySkew.wdl` exposes these inputs:

| Input | Type | Default | Description |
| --- | --- | --- | --- |
| `AnnotationData` | `File` | required | Annotated variant table. Plain TSV and gzip-compressed TSV inputs are supported. |
| `OutputFile` | `String` | required | Name for the aggregated output file, usually ending in `.tsv.gz`. |
| `VariantsPerShard` | `Int` | required | Number of input rows per shard. |
| `PipThreshold` | `Float` | `0.9` | Variants with `pip >= PipThreshold` are included. |
| `AdmixedSubpops` | `String` | `"oth"` | Comma-separated GVS subpopulation labels to remove for the no-admixed skew calculation. |

Example Cromwell-style inputs:

```json
{
  "ComputeAncestrySkew.AnnotationData": "eQTL_susie.COMB.annotated.cleaned.tsv.gz",
  "ComputeAncestrySkew.OutputFile": "AncestrySkew.tsv.gz",
  "ComputeAncestrySkew.VariantsPerShard": 50000,
  "ComputeAncestrySkew.PipThreshold": 0.9,
  "ComputeAncestrySkew.AdmixedSubpops": "oth"
}
```

### Required annotation columns

The R script requires:

| Column pattern | Description |
| --- | --- |
| `variant` | Variant identifier. |
| `pip` | Posterior inclusion probability used for filtering. |
| `gvs_all_ac`, `gvs_all_an` | Cohort-level alternate allele count and allele number. |
| `gvs_<subpop>_af` | Subpopulation allele frequency columns, for example `gvs_afr_af`, `gvs_amr_af`, `gvs_eur_af`. |
| `gvs_<subpop>_ac`, `gvs_<subpop>_an` | Matching subpopulation AC/AN columns for every `gvs_<subpop>_af` column. These are needed because max subpopulation is now recalculated from MAF and exact skew requires the matching AC/AN for the selected subpopulation. |

The script converts each `gvs_<subpop>_af` value to MAF with `min(AF, 1 - AF)` before choosing `gvs_max_subpop`. The existing input `gvs_max_subpop` and `gvs_max_af` columns are not trusted because they may have been computed from allele frequency rather than minor allele frequency.

## Outputs

The aggregated output is a gzip-compressed TSV with one row per variant passing the PIP threshold.

Key columns:

| Column | Description |
| --- | --- |
| `gvs_max_subpop` | MAF-derived max GVS subpopulation across all available subpopulations. |
| `gvs_max_maf` | MAF for `gvs_max_subpop`. |
| `gvs_max_ac`, `gvs_max_an` | AC/AN for `gvs_max_subpop`. |
| `gvs_background_ac`, `gvs_background_an` | `gvs_all_ac/an - gvs_max_ac/an`. |
| `gvs_odds_ratio`, `gvs_p_value` | Fisher test result for the standard ancestry skew calculation. |
| `gvs_no_admixed_all_ac`, `gvs_no_admixed_all_an` | Cohort AC/AN after subtracting the requested admixed subpopulation AC/AN. |
| `gvs_no_admixed_max_subpop` | MAF-derived max subpopulation after excluding the requested admixed subpopulations from max-subpopulation selection. |
| `gvs_no_admixed_background_ac`, `gvs_no_admixed_background_an` | No-admixed cohort AC/AN minus the no-admixed max subpopulation AC/AN. |
| `gvs_no_admixed_odds_ratio`, `gvs_no_admixed_p_value` | Fisher test result after admixed samples are removed from the background. |
| `gvs_no_admixed_removed_<subpop>_ac/an` | AC/AN subtracted for each requested admixed subpopulation. |

## Running The R Script Directly

```bash
Rscript scripts/ComputeAncestrySkew.R \
  --AnnotationData annotations.tsv.gz \
  --OutputPrefix output/annotations \
  --PipThreshold 0.9 \
  --AdmixedSubpops oth
```

This writes `output/annotations.AncestrySkew.tsv.gz`.
