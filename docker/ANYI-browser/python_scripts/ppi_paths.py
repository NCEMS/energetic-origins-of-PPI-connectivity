# ppi_paths.py (or include in ppi_io.py)
from __future__ import annotations

import os
from pathlib import Path


def resolve_data_dir(
    env_var: str = "PPI_DATA_DIR",
    repo_relative_default: str = "data",
) -> Path:
    """
    Resolve the directory containing nodes/edges data.

    Priority:
      1) $PPI_DATA_DIR if set
      2) ./data relative to the current working directory (typical repo root)
      3) raise a helpful error

    Returns a Path that exists.
    """
    # 1) environment variable (Docker-friendly)
    env_value = os.environ.get(env_var)
    if env_value:
        p = Path(env_value).expanduser().resolve()
        if p.exists():
            return p
        raise FileNotFoundError(
            f"{env_var} is set to '{env_value}', but that path does not exist."
        )

    # 2) repo-relative default (dev-friendly)
    p = (Path.cwd() / repo_relative_default).resolve()
    if p.exists():
        return p

    # 3) fail loudly with instructions
    raise FileNotFoundError(
        "Could not resolve data directory.\n"
        f"- Set {env_var} to the folder containing nodes.pkl and edges.csv, or\n"
        f"- Create a '{repo_relative_default}/' directory in your working directory.\n"
        f"Current working directory: {Path.cwd()}"
    )
