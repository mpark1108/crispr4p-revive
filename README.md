# CRISPR4P Revive

CRISPR4P Revive is a Python 3 primer- and donor-design tool for *Schizosaccharomyces pombe*, based on CRISPR4P. It supports sgRNA selection, cut-site annotation, deletion design, SpEDIT/pLSB cloning oligos, stop-cassette disruption, wild-type restoration, and PCR validation.

## Features

- Search by gene name or PomBase gene ID, chromosome interval, or oligo sequence
- Identify candidate sgRNAs and report genomic similarity counts
- Report PAM coordinates, Cas9 cut position, genomic region, CDS position, overlapping genes, and gene viability
- Generate CRISPR4P deletion-design primers
- Generate 52-nt SpEDIT/pLSB Golden Gate cloning oligos
- Generate disruption donors using a packaged set of 23-nt stop-cassette candidates with optional diagnostic sites (PacI or SwaI).
- Generate wild-type restoration donors and rescue-guide cloning oligos
- Design edit-spanning and optional junction-checking PCR primers

## Requirements

- Python 3
- primer3-py 2.3.0
- Gunicorn (Linux hosting)

The `primer3-py` package is imported in Python as `primer3`.

## Installation

```bash
git clone https://github.com/mpark1108/crispr4p-revive.git
cd crispr4p-revive
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
```

## Running the web application

```bash
python webapp.py
```

Then open [http://localhost:8080](http://localhost:8080).

### Hosting

On Linux, run the application with Gunicorn from the project directory:

```bash
.venv/bin/python -m gunicorn -c gunicorn.conf.py webapp:application
```

The supplied configuration uses four synchronous workers on `127.0.0.1:8080`.

## Reference data

The repository includes the *S. pombe* ASM294v2.26 reference genome and packaged PomBase gene-name, annotation, and gene-viability snapshots. Additional info in [`data/POMBASE_DATA.md`](data/POMBASE_DATA.md).

## Testing

```bash
python -m unittest discover -s tests -v
```

## License

CRISPR4P Revive is distributed under the MIT License. See [LICENSE](LICENSE).  

This project is derived from the original [Bähler Lab CRISPR4P](https://github.com/Bahler-Lab/crispr4p) software developed by Rodríguez-López, Fernández-Sánchez, and colleagues. The original CRISPR4P publication identifies the software as MIT-Licensed.

Rodríguez-López M, Cotobal C, Fernández-Sánchez O, et al.
“A CRISPR/Cas9-based method and primer design tool for seamless genome
editing in fission yeast.” *Wellcome Open Research*. 2017;1:19.

- [Journal article](https://wellcomeopenresearch.org/articles/1-19)
- [DOI: 10.12688/wellcomeopenres.10038.3](https://doi.org/10.12688/wellcomeopenres.10038.3)
- [![DOI](https://zenodo.org/badge/45244871.svg)](https://zenodo.org/badge/latestdoi/45244871)

## Related method

The SpEDIT/pLSB oligo workflow is based on:

Torres-Garcia S, Di Pompeo L, Eivers L, et al. “SpEDIT: A fast and efficient
CRISPR/Cas9 method for fission yeast.” *Wellcome Open Research*. 2020;5:274.

- [DOI: 10.12688/wellcomeopenres.16405.1](https://doi.org/10.12688/wellcomeopenres.16405.1)
