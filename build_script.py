"""Build rust backend.

This script is automatically run when the pyproject is installed.
"""

import csv
import io
import pathlib
import shutil
import subprocess


def download_gaia_data(output_file: pathlib.Path) -> None:
    """Download Gaia data and save it to the specified output file.

    Args:
        output_file: Path to save the downloaded Gaia data (csv).
    """
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = -1

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


def download_hipparcos_data(output_file: pathlib.Path) -> None:
    """Download Hipparcos data and save it to the specified output file.

    Args:
        output_file: Path to save the downloaded Hipparcos data (csv).
    """
    from astroquery.vizier import Vizier

    Vizier.ROW_LIMIT = -1

    print("Executing query on Hipparcos database...")
    query = Vizier.query_constraints(
        catalog=["I/311/hip2"],
        Hpmag="<7",
    )

    if not query:
        print("Job did not finish successfully.")
        return

    # Initialize a string buffer to write the CSV data to
    csv_output = io.StringIO()
    writer = csv.writer(csv_output)

    # Write the header row
    # Use column names that are common in Gaia and Hipparcos
    header = ["source_id", "ra", "dec", "parallax", "pmra", "pmdec", "phot_g_mean_mag"]
    writer.writerow(header)

    # Iterate over the rows of the Astropy table
    for row in query[0]:
        # Extract the required data points
        source_id = row["HIP"]
        ra = row["RArad"]
        dec = row["DErad"]
        parallax = row["Plx"]
        pmra = row["pmRA"]
        pmdec = row["pmDE"]
        # Use Vmag as the proxy for phot_g_mean_mag
        phot_g_mean_mag = row["Hpmag"]

        # Create a list for the new row and write it to the CSV writer
        new_row = [source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag]
        writer.writerow(new_row)

    # Get the complete CSV string from the buffer
    csv_string = csv_output.getvalue()

    # You can save this string to a file
    with output_file.open("w") as f:
        f.write(csv_string)

    print(f"Saved to: {output_file}")


def build_script() -> None:
    """Build rust backend and move shared library to correct folder."""
    cwd = pathlib.Path(__file__).parent.expanduser().absolute()

    gaia_file = cwd / "ruststartracker/gaia_data_j2016.csv"
    if not gaia_file.exists():
        download_gaia_data(gaia_file)

    hipparcos_file = cwd / "ruststartracker/hipparcos_data_j1991.25.csv"
    if not hipparcos_file.exists():
        download_hipparcos_data(hipparcos_file)

    subprocess.check_call(  # noqa: S603
        ["cargo", "build", "--release", "--features", "improc,gaia,hipparcos"],  # noqa: S607
        cwd=cwd,
        stdout=None,
        stderr=None,
    )
    shutil.copy(
        cwd / "target/release/libruststartracker.so", cwd / "ruststartracker/libruststartracker.so"
    )


if __name__ == "__main__":
    build_script()
