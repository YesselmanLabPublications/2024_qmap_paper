import logging
import sys

# logging #####################################################################

APP_LOGGER_NAME = "Q-DMS-TTR-PAPER"


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=False, file_name=None):
    """
    Set up the logger for the app

    Args:
        logger_name (str): The name of the logger (default is APP_LOGGER_NAME)
        is_debug (bool): Flag indicating whether to set the logger level to DEBUG (default is False)
        file_name (str): The name of the log file (default is None)

    Returns:
        logging.Logger: The configured logger instance
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)

    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    # pylint: disable=C0103
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)

    if file_name:
        # pylint: disable=C0103
        fh = logging.FileHandler(file_name)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


import logging


def get_logger(module_name):
    """
    Get the logger for the module.

    Args:
        module_name (str): The name of the module.

    Returns:
        logging.Logger: The logger object for the module.
    """
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)
