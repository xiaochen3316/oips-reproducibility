from pathlib import Path

import pandas as pd
import pytest


def load_io_module():
    try:
        from oips_repro import io
    except ImportError as exc:
        pytest.fail(f"deterministic I/O module is missing: {exc}")
    return io


def test_serialization_constants_are_frozen():
    io = load_io_module()
    assert io.ENCODING == "utf-8"
    assert io.LINE_TERMINATOR == "\n"
    assert io.FLOAT_FORMAT == "%.15g"
    assert io.NA_REP == ""
    assert io.BOOL_TRUE == "true"
    assert io.BOOL_FALSE == "false"


def test_write_stable_csv_sorts_and_serializes_deterministically(tmp_path: Path):
    io = load_io_module()
    frame = pd.DataFrame(
        {
            "key": ["b", "a"],
            "value": [1.2345678901234567, None],
            "flag": [True, False],
        }
    )
    output = tmp_path / "nested" / "stable.csv"

    io.write_stable_csv(
        frame,
        output,
        columns=["key", "value", "flag"],
        sort_by=["key"],
    )

    assert output.read_bytes() == (
        b"key,value,flag\n"
        b"a,,false\n"
        b"b,1.23456789012346,true\n"
    )
