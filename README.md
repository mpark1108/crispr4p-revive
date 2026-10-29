# CRISPR4P Revive

CRISPR4P Revive is a Python 3 modernization and extension of the original CRISPR4P primer-design tool for *Schizosaccharomyces pombe*. It preserves the original sgRNA and deletion-design workflow while adding cut-site annotation, SpEDIT/pLSB cloning oligos, stop-cassette disruption design, wild-type restoration design, and PCR validation primers.

## Features

- Search by gene name or PomBase gene ID, chromosome interval, or oligo sequence
- Identify candidate sgRNAs and report genomic similarity counts
- Report PAM coordinates, Cas9 cut position, genomic region, CDS position, overlapping genes, and gene viability
- Generate the original CRISPR4P deletion-design primers
- Generate 52-nt SpEDIT/pLSB Golden Gate cloning oligos
- Generate disruption donors using a packaged set of 23-nt stop-cassette candidates with optional diagnostic sites (PacI or SwaI).
- Generate wild-type restoration donors and rescue-guide cloning oligos
- Design edit-spanning and optional junction-checking PCR primers

## Requirements

- Python 3
- primer3-py

This version is currently tested with `primer3-py` 2.3.0. The `primer3-py` package is imported in Python as `primer3`.

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

### Apache CGI hosting

`webapp.py` runs as CGI when invoked by Apache. Running `python webapp.py`
directly starts the local preview server instead.

Configure an executable CGI wrapper named `webapp.py` to use the project's
virtual-environment Python, adjusting these paths to the installation:

```sh
#!/bin/sh
exec /path/to/crispr4p-revive/.venv/bin/python /path/to/crispr4p-revive/webapp.py
```

Serve `css/` as static files beside `webapp.py` in the public URL. Expose
only the CGI entry point and CSS, not the entire checkout.

The CGI account needs read access to the application, reference data, and
templates, and write access to `precomputed/`. This pickle cache must not
be writable by untrusted users.

The hosting administrator configures CGI execution, directory/SELinux
permissions, HTTPS, and service startup. See the
[Apache CGI guide](https://httpd.apache.org/docs/2.4/howto/cgi.html).

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
