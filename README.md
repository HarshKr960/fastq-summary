\---

title: "fastq-summary"
layout: page
---

# fastq-summary

[!\[CI](https://github.com/HarshKr960/fastq-summary/actions/workflows/ci.yml/badge.svg)](https://github.com/HarshKr960/fastq-summary/actions/workflows/ci.yml)
!\[Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)
!\[License](https://img.shields.io/badge/License-MIT-yellow.svg)

Dependency-free FASTQ summary statistics: read count, total bases, min/max/mean
length, GC%, mean Phred quality, and Q20/Q30 percentages. Plain-text or
gzip-compressed input, table or JSON output.

> Replace `YOUR-USERNAME` above with your actual GitHub username/org once pushed,
> or the badge will 404.

## Install

```bash
pip install .
# or, for development (adds pytest, flake8, coverage):
pip install -e ".\[dev]"
```

## Usage

```bash
fastq-summary reads.fastq.gz

fastq-summary reads\_R1.fastq.gz reads\_R2.fastq.gz --json
```

Table output:

```
file             num\_reads  total\_bases  min\_len  max\_len  mean\_len  gc\_pct  mean\_qual  q20\_pct  q30\_pct
---------------  ---------  -----------  -------  -------  --------  ------  ---------  -------  -------
reads.fastq.gz   1250000    45000000     36       36       36.0      41.2    34.8       97.1     91.4
```

## As a library

```python
from fastq\_summary import summarize

s = summarize("reads.fastq.gz")
print(s.gc\_content\_pct, s.mean\_quality)
print(s.to\_json())
```

## Testing

```bash
pip install -e ".\[dev]"
pytest --cov=fastq\_summary
```

Test fixtures live in `tests/data/` — a small hand-built FASTQ (`sample.fastq`,
plus a gzip copy) with statistics computed by hand in the test docstrings,
a length-mismatch file to exercise the error path, and an empty file.
No network access or real sequencing data needed, which is why this runs
cleanly in CI.

## CI/CD

`.github/workflows/ci.yml` runs on every push and pull request to `main`:

1. **Matrix test** across Python 3.9–3.12 — install, `flake8` lint, `pytest`
with coverage, upload the coverage report as a build artifact.
2. **CLI smoke test** — installs the built package and runs the actual
`fastq-summary` console script against the fixture file, asserting the
JSON output has the expected read count. This catches packaging/entry-point
breakage that unit tests alone wouldn't (e.g. a bad `pyproject.toml`
`\[project.scripts]` line).

## Known limitations

* Assumes Phred+33 quality encoding (correct for all modern Illumina data;
will silently misreport stats on legacy Phred+64 files).
* Loads per-file stats in memory as Python lists rather than streaming
aggregates — fine into the tens of millions of reads, but rewrite the
accumulation in `summarize()` to running sums if you need it on
hundred-million-read files under tight memory limits.

## License

MIT — see [LICENSE](LICENSE).

