"""Paths, session inventory and analysis-parameter loading.

Nothing in this module reads data; it only tells the rest of the package *where*
things are and *what* the analysis parameters are. Every tunable number lives in
``config/params.yaml`` and is reached through :func:`load_params` / :func:`param`.
"""

from __future__ import annotations

import os
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

#: The eight novel-maze sessions of CRCNS hc-11 (4 animals).
SESSIONS: tuple[str, ...] = (
    "Achilles_10252013",
    "Achilles_11012013",
    "Buddy_06272013",
    "Cicero_09012014",
    "Cicero_09102014",
    "Cicero_09172014",
    "Gatsby_08022013",
    "Gatsby_08282013",
)

#: Animal each session belongs to.
SESSION_ANIMAL: dict[str, str] = {s: s.split("_")[0] for s in SESSIONS}

#: Environment variable that overrides automatic data-directory discovery.
DATA_DIR_ENV_VAR = "HC11_DATA_DIR"

_SESSINFO_SUFFIX = "_sessInfo.mat"


def repo_root() -> Path:
    """Return the repository root (the directory containing ``pyproject.toml``)."""
    here = Path(__file__).resolve()
    for candidate in here.parents:
        if (candidate / "pyproject.toml").is_file():
            return candidate
    # Installed non-editably and detached from the source tree; fall back to cwd.
    return Path.cwd()


def _data_dir_candidates() -> list[Path]:
    root = repo_root()
    candidates: list[Path] = []
    env = os.environ.get(DATA_DIR_ENV_VAR)
    if env:
        candidates.append(Path(env).expanduser())
    candidates.append(root / "data" / "raw" / "NoveltySessInfoMatFiles")
    # Local convenience: data kept beside the repo rather than inside it.
    candidates.append(root.parent / "hc-11_Data" / "NoveltySessInfoMatFiles")
    return candidates


def data_dir(required: bool = True) -> Path:
    """Locate the directory holding the ``*_sessInfo.mat`` files.

    Resolution order: ``$HC11_DATA_DIR`` → ``<repo>/data/raw/NoveltySessInfoMatFiles``
    → ``<repo>/../hc-11_Data/NoveltySessInfoMatFiles``. The first candidate that
    exists *and* contains at least one ``*_sessInfo.mat`` file wins.

    Parameters
    ----------
    required
        If True (default), raise :class:`FileNotFoundError` when no candidate
        qualifies. If False, return the first candidate path regardless, so that
        callers can report a sensible location in an error message of their own.
    """
    candidates = _data_dir_candidates()
    for candidate in candidates:
        if candidate.is_dir() and any(candidate.glob(f"*{_SESSINFO_SUFFIX}")):
            return candidate
    if not required:
        return candidates[0]
    listed = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(
        "No hc-11 session files found. Looked in:\n  "
        f"{listed}\n"
        f"Set {DATA_DIR_ENV_VAR} to the directory containing the *{_SESSINFO_SUFFIX} "
        "files, or see the 'Getting the data' section of README.md."
    )


def session_path(session: str) -> Path:
    """Return the path to one session's ``sessInfo.mat``.

    Accepts either a bare session name (``"Achilles_10252013"``) or a full file
    name (``"Achilles_10252013_sessInfo.mat"``).
    """
    name = session[: -len(_SESSINFO_SUFFIX)] if session.endswith(_SESSINFO_SUFFIX) else session
    path = data_dir() / f"{name}{_SESSINFO_SUFFIX}"
    if not path.is_file():
        raise FileNotFoundError(f"Session file not found: {path}")
    return path


def available_sessions() -> list[str]:
    """Session names actually present on disk, in canonical order."""
    try:
        directory = data_dir()
    except FileNotFoundError:
        return []
    present = {p.name[: -len(_SESSINFO_SUFFIX)] for p in directory.glob(f"*{_SESSINFO_SUFFIX}")}
    known = [s for s in SESSIONS if s in present]
    extra = sorted(present - set(SESSIONS))
    return known + extra


def config_path(name: str = "params.yaml") -> Path:
    """Path to a file in the repository's ``config/`` directory."""
    return repo_root() / "config" / name


def results_dir(session: str | None = None, create: bool = True) -> Path:
    """Directory for generated figures and tables, optionally per session."""
    path = repo_root() / "results"
    if session is not None:
        path = path / session
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


@lru_cache(maxsize=8)
def _load_params_cached(path_str: str) -> dict[str, Any]:
    with open(path_str, encoding="utf-8") as fh:
        loaded = yaml.safe_load(fh)
    if not isinstance(loaded, dict):
        raise ValueError(f"{path_str} did not parse to a mapping.")
    return loaded


def load_params(path: str | Path | None = None) -> dict[str, Any]:
    """Load ``config/params.yaml`` (or another YAML file) as a nested dict.

    Results are cached per path; the returned dict is shared, so treat it as
    read-only.
    """
    resolved = Path(path) if path is not None else config_path()
    return _load_params_cached(str(resolved.resolve()))


def param(dotted_key: str, params: dict[str, Any] | None = None) -> Any:
    """Fetch a nested parameter by dotted path, e.g. ``"place_fields.bin_size_cm"``.

    Raises :class:`KeyError` naming the full path if any segment is missing, so a
    typo in a parameter name fails loudly instead of silently returning a default.
    """
    node: Any = load_params() if params is None else params
    for i, key in enumerate(dotted_key.split(".")):
        if not isinstance(node, dict) or key not in node:
            reached = ".".join(dotted_key.split(".")[:i]) or "<root>"
            raise KeyError(
                f"No parameter {dotted_key!r} in config (missing {key!r} under {reached})"
            )
        node = node[key]
    return node


def git_commit(short: bool = True) -> str:
    """Current git commit of the repository, or ``"unknown"`` outside a checkout.

    Appends ``"-dirty"`` when the working tree has uncommitted changes, so that a
    figure's sidecar JSON never claims a clean provenance it does not have.
    """
    try:
        args = ["git", "rev-parse", "--short" if short else "--verify", "HEAD"]
        commit = subprocess.check_output(
            args, cwd=repo_root(), stderr=subprocess.DEVNULL, text=True
        ).strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=repo_root(), stderr=subprocess.DEVNULL, text=True
        ).strip()
        return f"{commit}-dirty" if dirty else commit
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"


def _main() -> None:
    print(f"repo root    : {repo_root()}")
    print(f"git commit   : {git_commit()}")
    print(f"config       : {config_path()}  (exists={config_path().is_file()})")
    try:
        directory = data_dir()
        print(f"data dir     : {directory}")
    except FileNotFoundError as exc:
        print("data dir     : NOT FOUND")
        print(str(exc))
        return
    found = available_sessions()
    print(f"sessions     : {len(found)}/{len(SESSIONS)} found")
    for name in SESSIONS:
        mark = "ok     " if name in found else "MISSING"
        size = ""
        if name in found:
            size = f"{session_path(name).stat().st_size / 1e6:8.1f} MB"
        print(f"  [{mark}] {name:<20}{size}")
    for name in found:
        if name not in SESSIONS:
            print(f"  [unknown] {name}")


if __name__ == "__main__":
    _main()
