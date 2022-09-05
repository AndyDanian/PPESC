import os, sys

class HiddenPrints:
    """
    These class activate or desactivate the print function

    Based:
        https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
    Date:
        03/09/2022
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

def drv_hiddenprints(hidden_prints_other_object: HiddenPrints = None, verbose: int = -1):
    """
    Execute HiddenPrints object
    Args:
        hidden_prints_other_object (HiddenPrints): Oject that hidden or not prints
        verbose (int): one integer [-oo,0] off prints, one integer [1, oo] on print
    """

    if verbose <= 0:
        hidden_prints_other_object.__enter__()
    else:
        hidden_prints_other_object.__exit__(True,True,True)