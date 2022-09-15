from datetime import datetime
from pathlib import Path

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
    ##################################################################
    # PROPERTIES
    ##################################################################
    def scratch(self) -> Path:
        return self._scratch.parent
    def job_path(self) -> Path:
        return self._scratch
    ##################################################################
    # METHODS
    ##################################################################
    def remove_job_folder(self) -> None:
        """
        Remove job folder into scratch path
        """
        if self._scratch.exists() and self._scratch.is_dir:
            self._scratch.rmdir()
if __name__ == "__main__":
    s = scratch("/home1/scratch")
    print(s.scratch())
    print(s.job_path())
    s.remove_job_folder()
