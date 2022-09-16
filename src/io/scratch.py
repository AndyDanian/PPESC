from datetime import datetime
from pathlib import Path
import shutil
class scratch():
    # Constructor
    def __init__(self, scratch: Path or str = None, job_folder: str = None) -> None:
        """
        Constructor of scratch object

        Args:
            scrath (Path): Scrath's Path
        """

        if job_folder is None:
            now_date = datetime.now()
            str_date: str = now_date.strftime("%d%m%y%H%M%S")
            job_folder: str = str_date

        if scratch is None:
            if Path("/tmp").exists():
                tmp_scratch: Path = Path("/tmp/scratch").parent.mkdir(parents=True, exist_ok=True)
                self._scratch: Path = Path("tmp_scratch/"+(job_folder))
                self._scratch.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError("/tmp folder not exits.")
        else:
            if Path(scratch).exists():
                self._scratch: Path = Path(scratch + "/" + job_folder)
                if self._scratch.exists():
                    print("***WARNING\n\n{self._scratch} already exit, then possiblely the file overwrite")
                else:
                    self._scratch.mkdir(parents=True, exist_ok=True) 
            else:
                raise FileNotFoundError(f"***Error \n\n{Path(scratch)} folder not exists, please create before calculation.")
        
        self._output_name = None
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
    def output_name(self) -> str:
        return self._output_name
    @output_name.setter
    def output_name(self, name: str = None):
        if (self._scratch / (name)).exists():
            (self._scratch / (name)).rename((self._scratch / (name + ".0")))
        self._output_name = self._scratch / (name)

    ##################################################################
    # METHODS
    ##################################################################
    def write_header_output(self):
        """
        Write the header information for output file
        """
        with open(self._output_name, "a") as f:
            f.write(("*"*80)+"\n")
            f.write("PROGRAM NAME".center(80)+"\n")
            f.write("version 0.0".center(80)+"\n")
            f.write("2022".center(80)+"\n")
            f.write("Authors:".ljust(80)+"\n")
            f.write("         Andy Zapata".ljust(80)+"\n")
            f.write("*"*80+"\n")
            f.write("\n")

    def write_time(self, f: object = None, name: str = None, delta_time: float = None,
                header: bool = True, tailer: bool = True):
        """"
        Print time neccesary for calculations

        Args:
            name (str): Name of calculate
            delta_time (float): time in seconds
        """
        if len(name) > 60:
            count = 0
            words = ''
            for s in name.split():
                if count > 0:
                        count += 1
                count += len(s)
                if count <= 60:
                        words += s + ' '
                else:
                        count = 0
                        words += "\n" + s + ' '
            name = words

        if header:
            f.write("t"*20+"hours:minutes:seconds"+"t"*20)
        if delta_time <= 60:
            f.write(f"{name}, Time: 0:0:{delta_time:.3f}".ljust(62))
        elif delta_time > 60 and delta_time <= 3600:
            minutes = int(delta_time/60)
            seconds = delta_time%60
            f.write(f"{name}, Time: 0:{minutes}:{seconds:.3f}".ljust(62))
        else:
            hours = int(delta_time/3600)
            minutes = int(delta_time%3600/60)
            seconds = delta_time%3600%60
            f.write(f"{name}, Time: {hours}:{minutes}:{seconds:.3f}".ljust(62))
        if tailer:
            f.write("t"*63)

    def write_output(self, information: str = None, type: int = 0, delta_time: float = None,
                     header: bool = False, tailer: bool = False) -> list:
        """
        Save information for output file

        Args:
        ----
            information (str): Information to write in the output
            type (int): Type of information
                        0: standar information
                        1: time information
            delta_float (float): Delta of time for any process
            header (bool): Print header
            tailer (bool): Print tailer
        """


        if information is not None:      
            with open(self._output_name, "a") as f:
                if type == 0:
                    f.write(information+"\n")
                elif type == 1:
                    self.write_time(f, information, delta_time, header, tailer)

    def remove_job_folder(self) -> None:
        """
        Remove job folder into scratch path
        """
        if Path(self._output_name).exists():
            current_path: Path = Path().absolute() / (self._output_name.name)
            current_path.write_text(self._output_name.read_text())
        if self._scratch.exists() and self._scratch.is_dir:
            shutil.rmtree(self._scratch)
        
if __name__ == "__main__":
    s = scratch("/home1/scratch")
    print(s.scratch)
    print(s.job_path)
    s._output_name = "DATA.TXT"
    print("output name: ",s.output_name)
    s.remove_job_folder()
