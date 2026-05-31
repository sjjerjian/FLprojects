
from functools import wraps
import logging
from pathlib import Path
from time import time
from typing import Optional

def setup_loggers(
    log_file_path: Optional[Path | str] = None,
    log_level: Optional[int] = logging.INFO, # base level for all handlers
    print_to_console: bool = True,
    ):
    """Setup loggers for logging to console and file."""

    if log_file_path is not None:
        log_file_path = Path(log_file_path)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(log_level)

    formatter = logging.Formatter(
        '%(asctime)20s | %(levelname)7s | %(message)s', # %(module)20s:%(lineno)d - 
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # clear existing handlers
    for handler in logger.handlers:
        logger.removeHandler(handler)
        handler.close()

    if log_file_path is not None:
        file_handler = logging.FileHandler(log_file_path, mode='w')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    if print_to_console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)  
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    return logger


# decorator function to log time taken on function call
def time_fcn(func):
    @wraps(func)
    def wrap_func(*args, **kwargs):
        t1 = time.perf_counter()
        result = func(*args, **kwargs)
        t2 = time.perf_counter()
        logger.info(f'{func.__name__!r} executed in {(t2-t1):.4f}s')
        return result
    return wrap_func