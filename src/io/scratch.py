import sys
import os
from datetime import datetime
from pathlib import Path
import shutil
from typing import Union
from io import TextIOWrapper

import h5py
import numpy as np

from print_matrix import *

PROJECT_DIR = Path.cwd().parent
sys.path.append(os.fspath(PROJECT_DIR / "functions"))
from drv_time import drv_time

integral_symmetry: dict = {
    "overlap": "sym",
    "nucpot": "sym",
    "kinetic": "sym",
    "angmom": "antisym",
    "sd": "sym",
    "fc": "sym",
    "darwin": "sym",
    "massvelo": "sym",
    "nelfld": "sym",
    "diplen": "sym",
    "dipvel": "antisym",
    "pso": "antisym",
    "nstcgo": "sym",
    "dnske": "sym",
    "psoke": "square",
    "psooz": "square",
    "ozke": "antisym",
    "spinorbit": "antisym",
    "laplacian": "sym",
    "sofiel": "sym",
    "pnstcgop": "sym",
}


class scratch:
    # Constructor
    def __init__(self, scratch: Union[Path, str], job_folder: str) -> None:
        """
        Constructor of scratch object

        Args:
            scrath (Path): Scrath's Path
        """

        if not job_folder or job_folder == "":
            now_date = datetime.now()
            str_date: str = now_date.strftime("%d%m%y%H%M%S")
            job_folder = str_date

        if scratch is None:
            if Path("/tmp").exists():
                Path("/tmp/scratch").parent.mkdir(parents=True, exist_ok=True)
                job_folder = "/tmp/scratch/" + (job_folder)
                self._scratch: Path = Path(job_folder)
                self._scratch.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError("/tmp folder not exits.")
        else:
            if Path(scratch).exists():
                if isinstance(scratch, str):
                    scratch = Path(scratch)
                self._scratch = scratch / (job_folder)
                if self._scratch.exists():
                    print(
                        f"***WARNING\n\n{self._scratch} already exist, then possiblely the files will be overwrite"
                    )
                else:
                    self._scratch.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(
                    f"***Error \n\n{Path(scratch)} folder not exists, please create before calculation."
                )

        # H5 files
        h5_files: list = [
            "AO1BINT.H5",
            "AO2BINT.H5",
            "MO1BINT.H5",
            "EXCHCOUL.H5",
            "PRINPROP.H5",
        ]
        for i, name_file in enumerate(h5_files):
            if (self._scratch / (name_file)).exists():
                (self._scratch / (name_file)).unlink()
        self._hermite_ao1b_binary: Path = self._scratch / ("AO1BINT.H5")
        self._hermite_ao2b_binary: Path = self._scratch / ("AO2BINT.H5")
        self._hermite_mo1b_binary: Path = self._scratch / ("MO1BINT.H5")
        self._exchange_coulomb: Path = self._scratch / ("EXCHCOUL.H5")
        self._principal_propagator: Path = self._scratch / ("PRINPROP.H5")

        self._activate_write_output: bool = True

    ##################################################################
    # PROPERTIES
    ##################################################################
    @property
    def scratch(self) -> Path:
        return self._scratch.parent

    @property
    def job_path(self) -> Path:
        return self._scratch

    @property
    def output_path(self) -> Path:
        return self._output_path

    @output_path.setter
    def output_path(self, path: Path) -> None:
        if path.exists():
            path.rename(path.parent / (path.name + ".0"))
        self._output_path = path

    @property
    def activate_write_output(self) -> bool:
        return self._activate_write_output

    @activate_write_output.setter
    def activate_write_output(self, on_off: bool = True) -> None:
        self._activate_write_output = on_off

    ##################################################################
    # METHODS
    ##################################################################
    def write_header_output(self):
        """
        Write the header information for output file
        """
        with open(self._output_path, "a") as f:
            f.write(("*" * 80) + "\n")
            f.write("PROGRAM NAME".center(80) + "\n")
            f.write("version 0.0".center(80) + "\n")
            f.write("2022".center(80) + "\n")
            f.write("Authors:".ljust(80) + "\n")
            f.write("         Andy Zapata".ljust(80) + "\n")
            f.write("*" * 80 + "\n")
            f.write("\n")

    def write_box(self, f: TextIOWrapper, names: Union[str, list[str]], values: list):
        """
        Print a box with header and value

        Args:
            names (list): Names of the values to print
            values (list): Values of the values to print
        """

        l: int = len(names)
        lv: int = len(values)
        n: int = 1
        if l > 4:
            n += int(l / 4)

        if l > lv:
            raise ValueError("***ERROR\n\nThere are more header than values")
        elif lv % l != 0:
            raise ValueError("***ERROR\n\nValues size is not multiple of header size")

        values_lines: int = int(lv / l)

        split_up: str = ("┌" + "─" * 16 + "┐").center(18)
        split: str = ("├" + "─" * 16 + "┤").center(18)
        split_tailer: str = ("└" + "─" * 16 + "┘").center(18)

        for i in range(n):

            if i < n - 1:
                m: int = 4
            else:
                m = l - (n - 1) * 4

            header: str = ""
            split_ups: str = ""
            splits: str = ""
            split_tailers: str = ""
            for j in range(m):
                split_ups += split_up + " "
                splits += split + " "
                split_tailers += split_tailer + " "
                header += str("│" + "{}".format(names[j + i * 4]).center(16) + "│ ")

            f.write(split_ups.center(76) + "\n")
            f.write(header.center(76) + "\n")
            f.write(splits.center(76) + "\n")

            for k in range(values_lines):
                pvalues: str = ""
                for j in range(m):
                    if (
                        abs(values[j + i * 4 + k * l]) > 1.0e-2
                        and abs(values[j + i * 4 + k * l]) <= 9.9e6
                    ):
                        pvalues += str(
                            "│"
                            + "{:.3f}".format(values[j + i * 4 + k * l]).center(16)
                            + "│ "
                        )
                    else:
                        pvalues += str(
                            "│"
                            + "{:.5e}".format(values[j + i * 4 + k * l]).center(16)
                            + "│ "
                        )
                f.write(pvalues.center(76) + "\n")
            f.write(split_tailers.center(76) + "\n")

    def write_tensor(
        self,
        f: TextIOWrapper,
        names: Union[str, list[str]],
        values: list,
        isoani: bool = True,
        ani_axe: Union[str, int] = "z",
    ):
        """
        Print a matrix of 3x3 with or without iso/anisotropic

        Args:
            names (list): Names of the values to print
            values (list): Values of the values to print
            ani_axes (str or int): Axes to calculate the anisotropic value
            isoani (bool): Activate iso/anisotrpic print
        """
        sig_x: float = 1.0
        sig_y: float = 1.0
        sig_z: float = -1.0
        if ani_axe == 1 or ani_axe == 0 or ani_axe == "x":
            sig_x = -1.0
            sig_y = 1.0
            sig_z = 1.0
        elif ani_axe == 2 or ani_axe == "y":
            sig_x = 1.0
            sig_y = -1.0
            sig_z = 1.0

        l: int = len(names)
        lv: int = len(values)
        n: int = 1

        if l > lv:
            raise ValueError("***ERROR\n\nThere are more header than values")

        if l > 3:
            n += int(l / 3)

        split_up: str = ("┌" + "─" * 31 + "┐").center(32)
        split: str = ("├" + "─" * 31 + "┤").center(32)
        split_tailer: str = ("└" + "─" * 31 + "┘").center(32)

        for i in range(n):

            if i < n - 1:
                m: int = 3
            else:
                m = l - (n - 1) * 3

            header: str = ""
            split_ups: str = ""
            splits: str = ""
            split_tailers: str = ""

            for j in range(m):
                split_ups += split_up + " "
                splits += split + " "
                split_tailers += split_tailer + " "
                header += str("│" + "{}".format(names[j + i * 3]).center(31) + "│ ")
            f.write(split_ups.center(101) + "\n")
            f.write(header.center(101) + "\n")
            f.write(splits.center(101) + "\n")

            for k in range(3):
                pvalues: str = ""
                for j in range(m):
                    pvalues += "│ "
                    for p in range(3):
                        pvalues += str(
                            "{:.2e}".format(values[j + i * 3][k * 3 + p]).center(10)
                        )
                    pvalues += "│ "
                f.write(pvalues.center(101) + "\n")
            if isoani:
                f.write(splits.center(101) + "\n")
                iso_ani: str = ""
                for j in range(m):
                    iso_ani += str(
                        "│ISO: "
                        + "{:.3e}".format(
                            (
                                values[j + i * 3][0]
                                + values[j + i * 3][4]
                                + values[j + i * 3][8]
                            )
                            / 3.0
                        ).center(10)
                    )
                    iso_ani += str(
                        " ANI: "
                        + "{:.3e}".format(
                            (
                                sig_x * values[j + i * 3][0]
                                + sig_y * values[j + i * 3][4]
                                + sig_z * values[j + i * 3][8]
                            )
                            / 3.0
                        ).center(10)
                        + "│ "
                    )
                f.write(iso_ani.center(101) + "\n")
            f.write(split_tailers.center(101) + "\n")

    def write_title(self, f: TextIOWrapper, name: str = "", title_type: int = 0):
        """
        Print titles

        Args:
            name (str): Title name
        """
        if title_type == 0:
            title: str = "*** " + name.upper().center(70) + " ***"
            if len(name.upper()) > 70:
                print("***WARNING\n\n Title more large, please until 70 strings")
            f.write("*" * len(title) + "\n")
            f.write(title + "\n")
            f.write("*" * len(title) + "\n")
        elif title_type == 1:
            f.write("\n")
            f.write(name.center(70) + "\n")
            f.write(("-" * 40).center(70) + "\n")
        else:
            f.write(("=" * 40).center(80) + "\n")
            f.write(f"{name}".center(80) + "\n")
            f.write(("=" * 40).center(80) + "\n")

    def write_time(self, f: TextIOWrapper, drv_time: drv_time):
        """ "
        Print time neccesary for calculations

        Args:
            drv_time (object:drv_time): Driver the time process
        """

        count0 = 0
        for name, delta_time in zip(drv_time._name, drv_time._delta_time):
            header = False
            tailer = False
            if count0 == 0:
                header = True
            if count0 == len(drv_time._name) - 1:
                tailer = True
            count0 += 1

            if len(name) > 60:
                count = 0
                words = ""
                for s in name.split():
                    if count > 0:
                        count += 1
                    count += len(s)
                    if count <= 60:
                        words += s + " "
                    else:
                        count = 0
                        words += "\n" + s + " "
                name = words

            if header:
                f.write("t" * 20 + "hours:minutes:seconds" + "t" * 20 + "\n")
            if delta_time <= 60:
                f.write(f"{name}, Time: 0:0:{delta_time:.3f}".ljust(62) + "\n")
            elif delta_time > 60 and delta_time <= 3600:
                minutes = int(delta_time / 60)
                seconds = delta_time % 60
                f.write(f"{name}, Time: 0:{minutes}:{seconds:.3f}".ljust(62) + "\n")
            else:
                hours = int(delta_time / 3600)
                minutes = int(delta_time % 3600 / 60)
                seconds = delta_time % 3600 % 60
                f.write(
                    f"{name}, Time: {hours}:{minutes}:{seconds:.3f}".ljust(62) + "\n"
                )
            if tailer:
                f.write("t" * 63 + "\n")

    def write_size_file(
        self, f: TextIOWrapper, name_file: str = "", size_file: float = 0.0
    ):
        """
        Write size file into output

        Args:
        ----
            f (object): Object of output file
            size_file (float): Size of output file in bytes
        """
        unit: str = "Bytes"
        if size_file > 1024.0:
            size_file = size_file / (1024.0)
            unit = "KB"
        if size_file > 1024.0 * 1024.0:
            size_file = size_file / (1024.0 * 1024.0)
            unit = "MB"
        elif size_file > 1024.0 * 1024.0 * 1024.0:
            size_file = size_file / (1024.0 * 1024.0 * 1024.0)
            unit = "GB"

        f.write(f"{name_file}, size: {size_file:.3f} {unit} \n")

    def write_ao1bin_hermite(
        self,
        f: TextIOWrapper,
        array: dict[str, np.ndarray],
        direct: bool = False,
    ):
        """
        Print one body integrals into output

        Args:
        ----
            f (object): Object of output file
        """
        if not direct:
            with h5py.File(self._hermite_ao1b_binary, "r") as h:
                for name in list(h.keys()):
                    self.write_title(f, name, 1)
                    print_triangle_matrix(
                        f=f,
                        integral=h[name],
                        matriz_sym=integral_symmetry[name.split()[0]],
                    )
        else:
            for name, values in array.items():
                self.write_title(f, name, 1)
                if not name in integral_symmetry.keys():
                    symmetry: str = "square"
                else:
                    symmetry = integral_symmetry[name.split()[0].lower()]
                print_triangle_matrix(f=f, integral=values, matriz_sym=symmetry)

    def write_ao2bin_hermite(self, f: TextIOWrapper):
        """
        Print one body integrals into output

        Args:
        ----
            f (object): Object of output file
        """
        f.write("***Warning print two--integrals take a lot time")
        with h5py.File(self._hermite_ao2b_binary, "r") as h:
            for name in list(h.keys()):
                if name == "e2pot":
                    self.write_title(f, "Two--Body Repulsion Integrals", 1)
                else:
                    self.write_title(f, name, 1)
                for i in range(h[name].shape[0]):
                    for j in range(i + 1):
                        for k in range(i + 1):
                            if k < i:
                                m: int = k + 1
                            else:
                                m = j + 1
                            for l in range(m):
                                if abs(h[name][i, j, k, l]) > 1.0e-6:
                                    if abs(h[name][i, j, k, l]) > 999.0:
                                        formate: str = "{:.6e}"
                                    else:
                                        formate = "{:.6f}"
                                    f.write(
                                        f"{i+1:4} {j+1:4} {k+1:4} {l+1:4}    "
                                        + formate.format(h[name][i, j, k, l]).center(16)
                                        + "\n"
                                    )

    def write_output(
        self,
        information: str = "",
        type: int = 0,
        # title information
        title_type: int = 0,
        # time information
        drv_time: drv_time = None,
        # Size file in bytes
        size_file: float = 0.0,
        # Integral 1B
        direct=False,
        dictionary: dict[str, np.ndarray] = {},
        # Box
        box_type: int = 0,
        values: list = [],
    ) -> None:
        """
        Save information for output file

        Args:
        ----
            information (str): Information to write in the output
            type (int): Type of information
                        0: standar information
                        1: Titles
                            0: Main Title
                            1: Subtitle
                            2: Results
                        2: time information
                        3: size file
                        4: results into box
                            0: One box
                            1: Tensor box
                        9: hermite one body matriz
                       10: hermite two body integrals
            drv_time (object:drv_time): Driver of the time process
        """

        if self._activate_write_output:
            with open(self._output_path, "a") as f:
                if type == 0:
                    f.write(information + "\n")
                elif type == 1 or title_type > 0:
                    self.write_title(f, information, title_type)
                elif type == 2:
                    self.write_time(f, drv_time)
                elif type == 3:
                    self.write_size_file(f, information, size_file)
                elif type == 4:
                    if box_type == 0:
                        self.write_box(
                            f, information, values
                        )  #! change information by names type [str]
                    elif box_type == 1:  #!because this produce mypy error
                        self.write_tensor(f, information, values)
                elif type == 9:
                    self.write_ao1bin_hermite(f, dictionary, direct)
                elif type == 10:
                    self.write_ao2bin_hermite(f)

    def binary(
        self,
        file: Path = Path(),
        io: str = "",
        # Write information
        dictionary: dict = {},
        # Read or delete information
        label: str = "",
    ):
        """
        Save hermite information in AOINT.H5 binary file

        Args:
        ----
            file (Path): Path of binary file
            dictionary (dict): Information to write into binary file
            io (str): Indicate read:r, write:a, or delete:d in binary file
                    a: write
                    r: read
                    f: search
        Return:
        ------
            np.ndarray
        """
        if io is None or io.lower() not in ["a", "r", "d", "f"]:
            raise ValueError(
                f"***ERROR\n\n\
                            argument io due be a to write or r to read, io {io}"
            )

        if io.lower() == "d":
            wr: str = "a"
        elif io.lower() == "f":
            wr = "r"
        else:
            wr = io.lower()

        with h5py.File(file, wr.lower()) as h:
            if io.lower() == "a":
                for name, value in dictionary.items():
                    h[name] = value
            elif io.lower() == "r":
                return h[label][:]
            elif io.lower() == "f":
                return label in h.keys()
            else:
                del h[label]

    def remove_job_folder(self) -> None:
        """
        Remove job folder into scratch path
        """
        if Path(self._output_path).exists():
            current_path: Path = Path().absolute() / (self._output_path.name)
            current_path.write_text(self._output_path.read_text())
        if self._scratch.exists() and self._scratch.is_dir:
            shutil.rmtree(self._scratch)


if __name__ == "__main__":
    s = scratch("/home1/scratch", "")
    print(s.scratch)
    print(s.job_path)
    s._output_path = Path("DATA.TXT")
    print("output name: ", s.output_path)
    s.remove_job_folder()
