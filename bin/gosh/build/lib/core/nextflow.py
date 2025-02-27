import subprocess
from datetime import datetime

class NextflowRunner:
    def __init__(self):
        self.cmd = 'nextflow'
        
    def get_timestamp(self):
        """Get current timestamp in YYYYMMDD_HHMMSS format"""
        return datetime.now().strftime('%Y%m%d_%H%M%S')
        
    def run(self, args):
        """Run nextflow command with given arguments"""
        cmd = [self.cmd, 'run'] + args
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Pipeline execution failed: {e}")
