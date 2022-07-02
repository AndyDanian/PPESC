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
    Print titles

    Args:
        name (str): Subtitle name
    """
    print(name.title().center(70))
    print(("-" * 40).center(70))

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

