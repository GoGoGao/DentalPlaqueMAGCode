# Supragingival plaque multi-omics analysis

Scripts and core inputs for the Normal vs PoorPlaque shotgun metagenomics, metabolomics, and MAG analyses .

## Layout

```
code/           Figure-linked analysis scripts and paths_config.py
shared_input/   Core this-study omics matrices
```

Figure-to-script map: `FIGURE_SCRIPT_MAP.tsv`.

## `shared_input.zip`

| File | Description |
|------|-------------|
| `sample_group.tsv` | Sample IDs and Normal / PoorPlaque labels |
| `metaphlan_species_abundance.xls` | MetaPhlAn4 species relative abundances |
| `humann_ko_cpm.tsv` | HUMAnN KO abundances (CPM) |
| `metab_abundance_AB.tsv` | Metabolite abundances (n = 20) |
| `metab_desc.tsv` | Metabolite annotations |
| `mag_bin_abundance.tsv` | MAG relative abundances |
| `mag_gtdb_summary.tsv` | GTDB-Tk taxonomy and quality |
| `mag_ko_presence.xls` | MAG × KO matrix |
| `mag_phylogeny_unrooted.nwk` | MAG phylogeny |

Intermediate tables are not bundled and can be regenerated from these inputs.

## `code/`

Scripts are named by manuscript figure (e.g. `Fig1ABC_*.py`, `Fig5B_*.R`, `Fig6A_*.py`).  
Paths resolve via `paths_config.py` to `../shared_input/`.

