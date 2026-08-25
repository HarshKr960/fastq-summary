from pathlib import Path

import pytest

from fastq_summary.summary import FastqFormatError, iter_records, summarize

DATA = Path(__file__).parent / "data"


def test_iter_records_count():
    records = list(iter_records(DATA / "sample.fastq"))
    assert len(records) == 4
    header, seq, qual = records[0]
    assert header == "@read1 length10 highqual"
    assert seq == "ACGTACGTAC"
    assert len(seq) == len(qual)


def test_summarize_hand_computed_stats():
    """
    Values below are hand-computed from tests/data/sample.fastq:
      read1: len 10, GC=5,  qual all Phred 40
      read2: len 10, GC=2,  qual all Phred 0
      read3: len 8,  GC=8,  qual all Phred 20
      read4: len 12, GC=0,  qual all Phred 40
    """
    s = summarize(DATA / "sample.fastq")

    assert s.num_reads == 4
    assert s.total_bases == 40
    assert s.min_length == 8
    assert s.max_length == 12
    assert s.mean_length == pytest.approx(10.0)
    assert s.gc_content_pct == pytest.approx(37.5)
    assert s.mean_quality == pytest.approx(26.0)
    assert s.q20_pct == pytest.approx(75.0)
    assert s.q30_pct == pytest.approx(55.0)


def test_summarize_gzip_matches_plaintext():
    plain = summarize(DATA / "sample.fastq")
    gz = summarize(DATA / "sample.fastq.gz")

    assert plain.num_reads == gz.num_reads
    assert plain.total_bases == gz.total_bases
    assert plain.gc_content_pct == gz.gc_content_pct
    assert plain.mean_quality == gz.mean_quality


def test_malformed_seq_qual_length_mismatch_raises():
    with pytest.raises(FastqFormatError, match="length mismatch"):
        list(iter_records(DATA / "malformed.fastq"))


def test_empty_file_raises_on_summarize():
    with pytest.raises(FastqFormatError, match="No records found"):
        summarize(DATA / "empty.fastq")


def test_missing_file_raises_file_not_found():
    with pytest.raises(FileNotFoundError):
        summarize(DATA / "does_not_exist.fastq")


def test_to_dict_and_to_json_roundtrip():
    import json

    s = summarize(DATA / "sample.fastq")
    d = s.to_dict()
    assert d["num_reads"] == 4

    parsed = json.loads(s.to_json())
    assert parsed == d
