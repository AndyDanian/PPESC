def title(name: str = None):
    """
    Print titles

    Args:
        name (str): Title name
    """
    print("*" * 80)
    print("*** ", name.upper().center(70), " ***")
    print("*" * 80)