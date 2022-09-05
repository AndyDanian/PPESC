def print_title(name: str = None):
    """
    Print titles

    Args:
        name (str): Title name
    """
    print("*" * 80)
    print("*** ", name.upper().center(70), " ***")
    print("*" * 80)

def print_subtitle(name: str = None):
    """
    Print subtitles format

    Args:
        name (str): Subtitle name
    """
    print()
    print(name.center(70))
    print(("-" * 40).center(70))

def print_ljust(name:str = None):
    """
    Print left justyifique format

    Args:
        name (str): name
    """
    print()
    print(name.ljust(70))
    print(("-" * 40).ljust(70))

def print_result(name: str = None, value: float = None):
    """
    Print titles

    Args:
        name (str): Calculation name
        value (float, int): Calculation value
    """
    print(('='*40).ljust(70))
    print(f'{name}: {value}'.ljust(70))
    print(('='*40).ljust(70))

def print_box(names: list = None, values: list = None):
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
        n += int(l/4)

    if l > lv:
        raise ValueError("***ERROR\n\nThere are more header than values")
    elif lv%l != 0:
        raise ValueError("***ERROR\n\nValues size is not multiple of header size")

    values_lines: int = int(lv/l)

    split_up: str     = ("┌" + "─"*16 + "┐").center(18)
    split: str        = ("├" + "─"*16 + "┤").center(18)
    split_tailer: str = ("└" + "─"*16 + "┘").center(18)

    for i in range(n):

        if i < n-1:
            m: int = 4
        else:
            m: int = l - (n-1)*4

        header: str = ""
        split_ups: str = ""
        splits: str = ""
        split_tailers: str = ""
        for j in range(m):
            split_ups += split_up + " "
            splits += split + " "
            split_tailers += split_tailer + " "
            header += str("│" + "{}".format(names[j +i*4]).center(16) + "│ ")

        print(split_ups.center(76))
        print(header.center(76))
        print(splits.center(76))

        for k in range(values_lines):
            pvalues: str = ""
            for j in range(m):
                if abs(values[j + i*4 + k*l]) > 1.0e-2 and abs(values[j + i*4 + k*l]) <= 9.9e6:
                    pvalues += str("│" + "{:.3f}".format(values[j + i*4 + k*l]).center(16) + "│ ")
                else:
                    pvalues += str("│" + "{:.5e}".format(values[j + i*4 + k*l]).center(16) + "│ ")
            print(pvalues.center(76))
        print(split_tailers.center(76))

def print_tensor(names: list = None, values: list = None, isoani: bool = True, ani_axe: str or int = "z"):
    """
    Print a matrix of 3x3 with or without iso/anisotropic

    Args:
        names (list): Names of the values to print
        values (list): Values of the values to print
        ani_axes (str or int): Axes to calculate the anisotropic value
        isoani (bool): Activate iso/anisotrpic print
    """
    sig_x: int = 1.0
    sig_y: int = 1.0
    sig_z: int = -1.0
    if ani_axe == 1 or ani_axe == 0 or ani_axe == "x":
        sig_x: int = -1.0
        sig_y: int = 1.0
        sig_z: int = 1.0
    elif ani_axe == 2 or ani_axe == "y":
        sig_x: int = 1.0
        sig_y: int = -1.0
        sig_z: int = 1.0

    l: int = len(names)
    lv: int = len(values)
    n: int = 1

    if l > lv:
        raise ValueError("***ERROR\n\nThere are more header than values")

    if l > 3:
        n += int(l/3)

    split_up: str     = ("┌" + "─"*31 + "┐").center(32)
    split: str        = ("├" + "─"*31 + "┤").center(32)
    split_tailer: str = ("└" + "─"*31 + "┘").center(32)

    for i in range(n):

        if i < n-1:
            m: int = 3
        else:
            m: int = l - (n-1)*3

        header: str = ""
        split_ups: str = ""
        splits: str = ""
        split_tailers: str = ""

        for j in range(m):
            split_ups += split_up + " "
            splits += split + " "
            split_tailers += split_tailer + " "
            header += str("│" + "{}".format(names[j +i*3]).center(31) + "│ ")
        print(split_ups.center(101))
        print(header.center(101))
        print(splits.center(101))

        for k in range(3):
            pvalues: str = ""
            for j in range(m):
                pvalues += "│ "
                for p in range(3):
                    pvalues += str("{:.2e}".format(values[j + i*3][k*3 + p]).center(10))
                pvalues += "│ "
            print(pvalues.center(101))
        if isoani:
            print(splits.center(101))
            iso_ani: str = ""
            for j in range(m):
                iso_ani += str("│ISO: " +
                        "{:.3e}".format((values[j + i*3][0]+values[j + i*3][4]+values[j + i*3][8])/3.0).center(10))
                iso_ani += str(" ANI: " +
                        "{:.3e}".format((sig_x*values[j + i*3][0]+
                                        sig_y*values[j + i*3][4]+
                                        sig_z*values[j + i*3][8])/3.0).center(10)
                        + "│ ")
            print(iso_ani.center(101))
        print(split_tailers.center(101))

def print_time(name: str = None, delta_time: float = None,
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
        print()
        print("t"*20,"hours:minutes:seconds","t"*20)
    if delta_time <= 60:
        print(f"{name} Time: 0:0:{delta_time:.3f}".ljust(62))
    elif delta_time > 60 and delta_time <= 3600:
        minutes = int(delta_time/60)
        seconds = delta_time%60
        print(f"{name} Time: 0:{minutes}:{seconds:.3f}".ljust(62))
    else:
        hours = int(delta_time/3600)
        minutes = int(delta_time%3600/60)
        seconds = delta_time%3600%60
        print(f"{name} Time: {hours}:{minutes}:{seconds:.3f}".ljust(62))
    if tailer:
        print("t"*63,"\n")

