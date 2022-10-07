import sys
import os
import shutil
from pathlib import Path

import pytest

# Addres when execute from fock
RESPONSE_PATH = Path.cwd()
PARENT_PATH = RESPONSE_PATH.parent
sys.path.append(os.fspath(PARENT_PATH / ("io")))
from scratch import *


@pytest.mark.parametrize(
    "integrals_name, expected_results",
    [
        (
            [
                "overlap",
                "nucpot",
                "kinetic",
                "angmom",
                "sd",
                "fc",
                "darwin",
                "massvelo",
                "nelfld",
                "diplen",
                "dipvel",
                "pso",
                "nstcgo",
                "dnske",
                "psoke",
                "psooz",
                "ozke",
                "spinorbit",
                "laplacian",
                "sofiel",
                "pnstcgop",
            ],
            [
                "sym",
                "sym",
                "sym",
                "antisym",
                "sym",
                "sym",
                "sym",
                "sym",
                "sym",
                "sym",
                "antisym",
                "antisym",
                "sym",
                "sym",
                "square",
                "square",
                "antisym",
                "antisym",
                "sym",
                "sym",
                "sym",
            ],
        ),
    ],
)
def test_integral_symmetry_dict(integrals_name, expected_results):
    for count, name in enumerate(integrals_name):
        assert integral_symmetry[name] == expected_results[count]


def test_scratch_class_init():
    """
    Testing init in scratch class
    """
    s = scratch(scratch=Path().cwd())
    #! WARNINGS WITH THIS, BECAUSE YOU CAN ERRASE TESTS
    shutil.rmtree(s._scratch)
    assert s._scratch.parent == Path().cwd()


def test_scratch_tmp_folder():
    """
    Testing: the mkdir /test/scratch
    """
    s = scratch()
    #! WARNINGS WITH THIS, SHUTIL CAN ERRASE FILES OR FOLDERS
    #! PLEASE, YOU ARE SURE THE PATH
    shutil.rmtree(s._scratch.parent)
    assert s._scratch.parent == Path("/tmp/scratch")


def test_scratch_twice_tmp_folder(capfd):
    """
    Testing: the mkdir /test/scratch twice
    """
    s = scratch()
    sc = scratch(scratch=s._scratch.parent, job_folder=s._scratch.name)
    # * py.test --fixtures for a list of builtin fixtures.
    # * capfd fixture
    out, err = capfd.readouterr()
    #! WARNINGS WITH THIS, SHUTIL CAN ERRASE FILES OR FOLDERS
    #! PLEASE, YOU ARE SURE THE PATH
    shutil.rmtree(sc._scratch.parent)
    print("out 65 linea")
    assert (
        out
        == f"***WARNING\n\n{s._scratch} already exist, then possiblely the files will be overwrite\n"
    )


def test_FileNotFoundError_scratch():
    """
    Testing: the mkdir /home/scratch1 twice
    """
    with pytest.raises(FileNotFoundError):
        scratch(scratch=Path("/home/++aa99111**"))
