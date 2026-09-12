# Completing the classification of maximum scattered linear sets in PG(1,q^5)

This repository contains the preprint PDF and reproducibility materials for
the classification of maximum scattered linear sets in PG(1,q^5).

## Paper

- [Preprint PDF, including appendices](paper/classification_maximum_scattered_linear_sets_PG1q5.pdf).

Archived preprint: [DOI 10.5281/zenodo.22232593](https://doi.org/10.5281/zenodo.22232593).

## Repository structure

| Location | Contents |
| --- | --- |
| [paper/](paper/) | Preprint and supplementary appendix PDFs. |
| [code/](code/README.md) | Mathematical verification programs, dependencies, and instructions. |
| [data/](data/README.md) | Required inputs and reference numerical results. |
| [CITATION.cff](CITATION.cff) | Machine-readable author and citation metadata. |
| [LICENSE](LICENSE) | License and copyright notice. |
| [manifest.json](manifest.json) | File inventory, sizes, and SHA-256 digests. |
| [SHA256SUMS.txt](SHA256SUMS.txt) | SHA-256 checksums for this public snapshot. |

## Reproduce the checks

Install the locked environment described in [code/README.md](code/README.md), then run from the repository root:

```text
python code/run_checks.py
```

For the locked environment, use `uv run --project code --locked python -B code/run_checks.py`.

See [verification instructions](code/README.md) for installation, optional
modes, expected results, and runtime. [Data notes](data/README.md) describe
the inputs and reference outputs.

## Scope of the checks

The programs check symbolic matrix and Pfaffian identities, finite-field calculations, and direct collision certificates. The full command enumerates admissible C3 and C4 parameters for every prime power q <= 25. These finite computations supplement the proof of the general classification.

The mathematical proofs are contained in the PDFs. This preprint has not
undergone journal peer review.

## Citation

Please cite the article using [CITATION.cff](CITATION.cff).

## License

**All rights reserved**

© 2026 Shunqi Lu. All rights reserved.

See [LICENSE](LICENSE).

## Author

Shunqi Lu, School of Artificial Intelligence, Capital University of Economics
and Business, Beijing 100070, China.

- ORCID: https://orcid.org/0009-0005-8251-2070
- Correspondence: 32023210246@cueb.edu.cn
