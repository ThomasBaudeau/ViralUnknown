from pathlib import Path
import subprocess

def check_dependencies():
    """Check that the required tools are available."""
    tools = ["snakemake", "conda"]
    missing = []
    for tool in tools:
        result = subprocess.run(
            ["which", tool],
            capture_output=True
        )
        if result.returncode != 0:
            missing.append(tool)
    return missing