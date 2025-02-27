import csv
import os

def check_if_tumor_only(samplesheet_path):
    """
    Check the 'status' column in the samplesheet to determine if the pipeline should run in tumor-only mode.
    Tumor-only mode is enabled if there are no rows with 'status' == '0' (i.e., no normal samples).

    Args:
        samplesheet_path (str): Path to the samplesheet CSV file.

    Returns:
        bool: True if tumor-only mode should be enabled, False otherwise.
    """
    has_normal = False
    if not os.path.exists(samplesheet_path):
        raise FileNotFoundError(f"Samplesheet not found at {samplesheet_path}")

    with open(samplesheet_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        if 'status' not in reader.fieldnames:
            raise ValueError("Samplesheet missing required 'status' column.")
        for row in reader:
            if row.get('status') == '0':
                print(f"Normal sample found: {row}")
                has_normal = True
                break
    return not has_normal  # Returns True if tumor-only (no normal samples found)
