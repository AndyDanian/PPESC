from string_informations import *


class drv_time:
    def __init__(self):
        """
        Object with the resposibility of the administration of the time print
        """

        self._name: list[str] = []
        self._delta_time: list[float] = []

    ##################################################################
    # ATTRIBUTES
    ##################################################################
    @property
    def name(self) -> list[str]:
        "Sections Name"
        return self._name

    @property
    def delta_time(self) -> list[float]:
        "Sections Time"
        return self._delta_time

    @property
    def save_number(self) -> int:
        "Sections Number"
        return len(self._name)

    @property
    def print_information(self) -> None:
        return print(
            [
                name + " time " + str(t) + " seconds "
                for name, t in zip(self._name, self._delta_time)
            ]
        )

    @property
    def reset(self):
        self._name: list[str] = []
        self._delta_time: list[float] = []

    ##################################################################
    # METHODS
    ##################################################################
    def add_name_delta_time(self, name: str, delta_time: float) -> None:
        """
        Save name and time delta in the respective variables

        Args:
        ----
            name (str): Name or phrase that descripte calculation
            delta_time (float): Calculation time
        """
        if name is None:
            raise ValueError("***ERROR\n\n\ Section name is neccesary.")
        if delta_time is None:
            raise ValueError("***ERROR\n\n\ Time delta is neccesary.")
        self._name.append(name)
        self._delta_time.append(delta_time)

    def printing(self) -> None:
        """
        Print the information
        """

        count = 0
        for name, dtime in zip(self._name, self._delta_time):
            header = False
            tailer = False
            if count == 0:
                header = True
            if count == self.save_number - 1:
                tailer = True
            print_time(name=name, delta_time=dtime, header=header, tailer=tailer)
            count += 1


if __name__ == "__main__":
    drv_t = drv_time()
    drv_t.add_name_delta_time(name="AA EE", delta_time=0.3)
    drv_t.add_name_delta_time(name="FF QQ", delta_time=1.0)
    drv_t.printing()
