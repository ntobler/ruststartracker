"""Build rust backend.

This script is automatically run when the pyproject is installed.
"""

import pathlib
import shutil
import subprocess


def rust_build() -> None:
    """Build rust backend and move shared library to correct folder."""
    cwd = pathlib.Path(__file__).parent.expanduser().absolute()
    subprocess.check_call(["cargo", "build", "--release"], cwd=cwd)  # noqa: S603, S607
    shutil.copy(
        cwd / "target/release/libruststartracker.so", cwd / "ruststartracker/libruststartracker.so"
    )
    subprocess.check_call(["cargo", "clean"], cwd=cwd)  # noqa: S603, S607


if __name__ == "__main__":
    rust_build()
