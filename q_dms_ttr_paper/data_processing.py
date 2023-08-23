import pandas as pd

from q_dms_ttr_paper.logger import setup_applevel_logger, get_logger


import pandas as pd
from pathlib import Path


def get_preprocessed_data(path, sets) -> pd.DataFrame:
    dfs = []
    for run_name in sets:
        full_path = Path(path) / run_name / "analysis" / "summary.json"
        if not full_path.exists():
            raise FileNotFoundError(f"File {full_path} does not exist")
        df = pd.read_json(full_path)
        dfs.append(df)
    return pd.concat(dfs)
