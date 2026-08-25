import json
from pathlib import Path

from fastq_summary.cli import main

DATA = Path(__file__).parent / "data"


def test_cli_table_output(capsys):
    rc = main([str(DATA / "sample.fastq")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "num_reads" in out
    assert "4" in out


def test_cli_json_output(capsys):
    rc = main([str(DATA / "sample.fastq"), "--json"])
    assert rc == 0
    out = capsys.readouterr().out
    parsed = json.loads(out)
    assert parsed[0]["num_reads"] == 4


def test_cli_multiple_files(capsys):
    rc = main([str(DATA / "sample.fastq"), str(DATA / "sample.fastq.gz")])
    assert rc == 0
    out = capsys.readouterr().out
    assert out.count("sample.fastq") >= 2


def test_cli_error_exit_code(capsys):
    rc = main([str(DATA / "empty.fastq")])
    assert rc == 1
    err = capsys.readouterr().err
    assert "error:" in err
