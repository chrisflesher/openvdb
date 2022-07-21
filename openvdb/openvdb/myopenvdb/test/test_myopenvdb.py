"""
Unit tests for the OpenVDB Python module

These are intended primarily to test the Python-to-C++ and
C++-to-Python bindings, not the OpenVDB library itself.
"""

import pytest

import myopenvdb as vdb


@pytest.mark.parametrize('level', [
    # 'debug',
    # 'info',
    'warn',
    # 'error',
    # 'fatal',
])
def test_set_logger_level(level):
    vdb.setLoggingLevel(level)
    assert vdb.getLoggingLevel() == level


def test_set_invalid_logger_level():
    with pytest.raises(ValueError):
        vdb.setLoggingLevel('invalid')


def test_set_program_name():
    vdb.setProgramName('name')
    vdb.setProgramName('name', True)
    vdb.setProgramName('name', False)
