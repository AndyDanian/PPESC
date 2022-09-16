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
                self._scratch: Path = Path(tmp_scratch + "/" + job_folder)
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

    def write_output(self, information: str = None) -> list:
        """
        Save information for output file

        Args:
        ----
            information (str): Information to write in the output
        """
        if information is None:
            raise ValueError("***ERROR\n\n\ Section name is neccesary.")
        
        with open(self._output_name, "a") as f:
            f.write(information)
        

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
