"""
Command-line interface: fastq-summary file1.fastq[.gz] [file2 ...]
"""

from __future__ import annotations

import argparse
import sys

from .summary import FastqFormatError, summarize


def _print_table(summaries) -> None:
    headers = [
        "file", "num_reads", "total_bases", "min_len", "max_len",
        "mean_len", "gc_pct", "mean_qual", "q20_pct", "q30_pct",
    ]
    rows = [
        [
            s.file, s.num_reads, s.total_bases, s.min_length, s.max_length,
            s.mean_length, s.gc_content_pct, s.mean_quality, s.q20_pct, s.q30_pct,
        ]
        for s in summaries
    ]
    widths = [
        max(len(str(h)), max((len(str(r[i])) for r in rows), default=0))
        for i, h in enumerate(headers)
    ]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print(fmt.format(*["-" * w for w in widths]))
    for r in rows:
        print(fmt.format(*r))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="fastq-summary",
        description="Summarize FASTQ read count, length, GC%, and quality.",
    )
    parser.add_argument("files", nargs="+", help="FASTQ file(s), .fastq/.fq, optionally .gz")
    parser.add_argument(
        "--json", action="store_true", help="Emit JSON instead of a table"
    )
    args = parser.parse_args(argv)

    summaries = []
    for f in args.files:
        try:
            summaries.append(summarize(f))
        except (FastqFormatError, FileNotFoundError) as exc:
            print(f"error: {f}: {exc}", file=sys.stderr)
            return 1

    if args.json:
        import json
        print(json.dumps([s.to_dict() for s in summaries], indent=2))
    else:
        _print_table(summaries)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
