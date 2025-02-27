import os
import click
from core.nextflow import NextflowRunner
from core.config import Config
from core.params_wizard import create_params_file

@click.group(name='run')
def run_cli():
    """Run pipeline commands"""
    pass

@run_cli.command()
@click.option('--pipeline-dir',
              default='/gpfs/commons/home/sdider/Projects/nf-casereports',
              help='Path to pipeline directory')
@click.option('--params-file',
              default='./params.json',
              help='Path to parameters JSON file')
@click.option('--profile',
              default='singularity',
              help='Pipeline profile(s) to use; comma-separated')
@click.option('--resume/--no-resume',
              default=True,
              help='Resume previous run if possible')

def pipeline(pipeline_dir, params_file, profile, resume):
    """Run the nextflow pipeline with specified parameters"""

    # Create hg19 directory (required for fragcounter)
    os.makedirs('hg19', exist_ok=True)

    # Check if params_file is provided and exists
    if not os.path.isfile(params_file):
        print(f"Parameters file '{params_file}' not found.")
        # Check if default params.json exists
        default_params = './params.json'
        if not os.path.isfile(default_params):
            print("No params.json file found. Launching wizard to create one...")
            # Call the wizard to create params.json
            create_params_file()
            params_file = default_params
        else:
            params_file = default_params
            print(f"Using default parameters file '{params_file}'.")

    # Initialize runner
    runner = NextflowRunner()

    # Build command arguments
    args = [
        pipeline_dir,
        '-params-file', params_file,
        '-profile', profile,
        '-with-report', f"report_{runner.get_timestamp()}.html",
        '-with-trace'
    ]

    if resume:
        args.append('-resume')

    # Execute pipeline
    runner.run(args)
