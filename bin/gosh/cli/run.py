import os
import click
from core.nextflow import NextflowRunner
from core.config import Config

@click.group(name='run')
def run_cli():
    """Run pipeline commands"""
    pass

@run_cli.command()
@click.option('--pipeline-dir',
              default='/gpfs/commons/home/sdider/Projects/nf-casereports',
              help='Path to pipeline directory')
@click.option('--nygc-config',
              default='/gpfs/commons/home/sdider/Projects/nf-casereports/tests/test_runs/nygc.config',
              help='Path to NYGC config file')
@click.option('--params-file',
              default='./params.json',
              help='Path to parameters JSON file')
@click.option('--profile',
              default='singularity,chr21_test',
              help='Pipeline profile to use')
@click.option('--resume/--no-resume',
              default=True,
              help='Resume previous run if possible')

def pipeline(pipeline_dir, nygc_config, params_file, profile, resume):
    """Run the nextflow pipeline with specified parameters"""

    # Set environment variables
    os.environ['NXF_SINGULARITY_LIBRARYDIR'] = '/gpfs/commons/groups/imielinski_lab/data/pipeline/container_images_cache'
    os.environ['NXF_SINGULARIY_CACHEDIR'] = os.path.expanduser('~/local_singularity_cache')

    # Create cache directory if it doesn't exist
    os.makedirs(os.environ['NXF_SINGULARIY_CACHEDIR'], exist_ok=True)

    # Create hg19 directory (required for fragcounter)
    os.makedirs('hg19', exist_ok=True)

    # Initialize runner
    runner = NextflowRunner()

    # Build command arguments
    args = [
        pipeline_dir,
        '-c', f"{pipeline_dir}/tests/nextflow_chr21_test.config",
        '-c', nygc_config,
        '-dump-channels',
        '-params-file', params_file,
        '-profile', profile,
        '-with-report', f"report_{runner.get_timestamp()}.html",
        '-with-trace'
    ]

    if resume:
        args.append('-resume')

    # Execute pipeline
    runner.run(args)
