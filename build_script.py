"""Build rust backend.

This script is automatically run when the pyproject is installed.
"""

import pathlib
import shutil
import subprocess


def download_gaia_data(output_file: pathlib.Path) -> None:
    """Download Gaia data and save it to the specified output file.

    Args:
        output_file: Path to save the downloaded Gaia data (csv).
    """
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = 10000

    query = """
    SELECT source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag
    FROM gaiadr3.gaia_source
    WHERE phot_g_mean_mag < 7
        AND parallax > 0
        AND astrometric_params_solved = 31
        AND visibility_periods_used >= 8
    """

    print("Executing query on Gaia database...")
    job = Gaia.launch_job_async(
        query,
        dump_to_file=True,
        output_format="csv",
        output_file=str(output_file.absolute()),
        verbose=True,
    )
    if not job.is_finished():
        print("Job did not finish successfully.")
        return

    print(f"Saved to: {output_file}")


def build_script() -> None:
    """Build rust backend and move shared library to correct folder."""
    cwd = pathlib.Path(__file__).parent.expanduser().absolute()

    gaia_file = cwd / "ruststartracker/gaia_data_j2016.csv"
    if not gaia_file.exists():
        download_gaia_data(gaia_file)

    subprocess.check_call(  # noqa: S603
        ["cargo", "build", "--release", "--features", "improc,gaia"],  # noqa: S607
        cwd=cwd,
        stdout=None,
        stderr=None,
    )
    shutil.copy(
        cwd / "target/release/libruststartracker.so", cwd / "ruststartracker/libruststartracker.so"
    )


if __name__ == "__main__":
    build_script()
