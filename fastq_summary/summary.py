"""
fastq_summary.summary
----------------------
Streaming FASTQ parser and per-file summary statistics.

Supports plain-text and gzip-compressed FASTQ (.fastq, .fq, .fastq.gz, .fq.gz).
No third-party dependencies — stdlib only, so CI has nothing extra to install.
"""

from __future__ import annotations

import gzip
import io
import json
import statistics
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterator, Tuple, Union

PathLike = Union[str, Path]

# Phred+33 is by far the dominant encoding on modern Illumina data
# (Phred+64 was retired with the GA IIx / early HiSeq era).
PHRED_OFFSET = 33


class FastqFormatError(ValueError):
    """Raised when a file does not look like well-formed FASTQ."""


def _open_maybe_gzip(path: PathLike) -> io.TextIOBase:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_records(path: PathLike) -> Iterator[Tuple[str, str, str]]:
    """
    Yield (header, sequence, quality) tuples for each FASTQ record.

    Raises FastqFormatError on truncated records or a header line that
    doesn't start with '@'.
    """
    with _open_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return  # clean EOF between records
            header = header.rstrip("\n")
            if not header.startswith("@"):
                raise FastqFormatError(
                    f"Expected header line starting with '@', got: {header!r}"
                )

            seq = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            qual = handle.readline().rstrip("\n")

            if not plus.startswith("+"):
                raise FastqFormatError(
                    f"Expected '+' separator line, got: {plus!r}"
                )
            if len(seq) != len(qual):
                raise FastqFormatError(
                    f"Sequence/quality length mismatch in record {header!r}: "
                    f"{len(seq)} bases vs {len(qual)} quality scores"
                )

            yield header, seq, qual


@dataclass
class FastqSummary:
    file: str
    num_reads: int
    total_bases: int
    min_length: int
    max_length: int
    mean_length: float
    gc_content_pct: float
    mean_quality: float
    q20_pct: float
    q30_pct: float

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)


def summarize(path: PathLike) -> FastqSummary:
    """Compute summary statistics for a single FASTQ file."""
    lengths = []
    gc_count = 0
    base_count = 0
    qual_scores = []
    q20 = 0
    q30 = 0

    for _, seq, qual in iter_records(path):
        lengths.append(len(seq))
        base_count += len(seq)
        gc_count += sum(1 for b in seq.upper() if b in "GC")
        for ch in qual:
            score = ord(ch) - PHRED_OFFSET
            qual_scores.append(score)
            if score >= 20:
                q20 += 1
            if score >= 30:
                q30 += 1

    num_reads = len(lengths)
    if num_reads == 0:
        raise FastqFormatError(f"No records found in {path}")

    total_quality_calls = len(qual_scores)

    return FastqSummary(
        file=str(path),
        num_reads=num_reads,
        total_bases=base_count,
        min_length=min(lengths),
        max_length=max(lengths),
        mean_length=round(statistics.mean(lengths), 2),
        gc_content_pct=round(100 * gc_count / base_count, 2) if base_count else 0.0,
        mean_quality=round(statistics.mean(qual_scores), 2) if qual_scores else 0.0,
        q20_pct=round(100 * q20 / total_quality_calls, 2) if total_quality_calls else 0.0,
        q30_pct=round(100 * q30 / total_quality_calls, 2) if total_quality_calls else 0.0,
    )
