import subprocess
from datetime import datetime

class NextflowRunner:
    def __init__(self):
        self.cmd = 'nextflow'

    def get_timestamp(self):
        """Get current timestamp in YYYYMMDD_HHMMSS format"""
        return datetime.now().strftime('%Y%m%d_%H%M%S')

    def run(self, command):
        """Run nextflow command with given command string"""
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Pipeline execution failed: {e}")
