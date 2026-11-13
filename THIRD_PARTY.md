# Upstream tools and data

| Component | Upstream | Use / release action |
|---|---|---|
| PlantCAD / PlantCAD2 | https://github.com/plantcad/plantcad | Backbone and methodological source; Apache-2.0; license retained in licenses/PlantCAD-LICENSE |
| GeneCAD | https://github.com/plantcad/genecad | v0.4.0, locally resolved commit 1a9baf84d61152c20a0cd910694316cd88d1d374; Apache-2.0; license retained in licenses/GeneCAD-LICENSE; see gene_annotation/UPSTREAM_VERSION.md |
| Evo2 | https://github.com/ArcInstitute/evo2 | Autoregressive benchmark; retain directional decoding and exact checkpoint identifiers |
| OrthoFinder | https://github.com/davidemms/OrthoFinder | Version 2.5.4, default settings; reference PMID 31727128 |
| Marchantia variants | https://marchantia.info/variants/data/Marpol_135_accessions/ | Natural SNP dataset, 135 accessions; PMID 39962240 |
| M. polymorpha accessibility | GEO GSE290177; PMID 41391150 | Locally processed/merged peaks still needed |
| A. thaliana accessibility | GEO GSE145742; PMID 32526034 | DNase-seq, not ATAC-seq; locally processed peaks still needed |

The 126-assembly orthogroup resource and the 135-accession variant dataset are distinct. Dataset citations and assembly provenance are recorded in manifests. Metadata listings do not grant redistribution rights to full genome files.

Unchanged upstream packages should be installed from fixed versions. For modified tools, distribute the matching patch/fork and retain upstream notices. BryoCAD includes scripts developed or adapted using PlantCAD/PlantCAD2 and GeneCAD functionality. Original notices are retained; GeneCAD modifications are documented with the pinned patch. The SNP scorer includes BryoCAD-specific input-handling fixes documented in variant_analysis/SCORING_REVIEW.md.

The root Apache-2.0 license covers BryoCAD code. Genome sequence excerpts and other scientific input data retain their source terms; model weights are governed by their model cards.

The unified environment builds [Mamba v2.2.4](https://github.com/state-spaces/mamba/tree/v2.2.4), [causal-conv1d v1.5.0.post8](https://github.com/Dao-AILab/causal-conv1d/tree/v1.5.0.post8), and FlashAttention 2.7.4.post1. The architecture-flag patch extends the author-provided dicos build procedure; exact native source commits, setup.py hashes and diffs are retained in build records. These dependencies retain their upstream licenses.
