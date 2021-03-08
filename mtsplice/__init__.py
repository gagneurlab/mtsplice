import os


def get_data_dir():
    """Returns the data directory
    """
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    this_path = os.path.dirname(os.path.abspath(filename))
    DATA = os.path.join(this_path, "../data/")
    if not os.path.exists(DATA):
        raise ValueError(DATA + " folder doesn't exist")
    return DATA

def get_paper_dir():
    """Returns the directory that stored models
    """
    import inspect
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    this_path = os.path.dirname(os.path.abspath(filename))
    PAPER = os.path.join(this_path, "../paper/")
    if not os.path.exists(PAPER):
        raise ValueError(PAPER + " folder doesn't exist")
    return PAPER


DATADIR = get_data_dir()
