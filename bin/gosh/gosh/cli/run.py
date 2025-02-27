import os
import click
import socket
import shutil
import subprocess
import sys
from ..core.nextflow import NextflowRunner
from ..core.config import Config
from ..core.params_wizard import create_params_file

@click.group(name='run')
def run_cli():
    """Run pipeline commands"""
    pass

def get_environment_defaults():
    nyu_defaults = {
        'pipeline-dir': "/gpfs/data/imielinskilab/projects/nf-casereports",
        'profile': "nyu",
        'nextflow_module': "nextflow/23.04.4"
    }

    mapping = {
        'cn-': nyu_defaults,
        'gn-': nyu_defaults,
        'bigpurple': nyu_defaults,
    }

    hostname = socket.gethostname()

    for prefix, defaults in mapping.items():
        if hostname.startswith(prefix):
            print(f"Detected environment: {hostname}")
            print(f"Using defaults for {hostname}")
            return defaults
    return {}

def load_required_modules(env_defaults):
    """Load required modules if commands are not available."""
    required_commands = ['nextflow', 'aws']
    modules_to_load = []
    load_modules_command = ""

    # Check for 'nextflow' command
    if shutil.which('nextflow') is None:
        nextflow_module = env_defaults.get('nextflow_module', 'nextflow')
        if not nextflow_module:
            nextflow_module = 'nextflow'
        modules_to_load.append(nextflow_module)
        print(f"'nextflow' command not found. Loading module '{nextflow_module}'.")
    else:
        print("'nextflow' command is already available.")

    # Check for 'aws' command
    if shutil.which('aws') is None:
        modules_to_load.append('aws-cli')
        print("'aws' command not found. Loading module 'aws-cli'.")
    else:
        print("'aws' command is already available.")

    # Load required modules using 'modulecmd'
    for module in modules_to_load:
        load_modules_command += f"module load {module} && "

    return load_modules_command

@run_cli.command()
@click.option('--pipeline-dir',
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

    # Retrieve environment defaults
    env_defaults = get_environment_defaults()

    # Set pipeline_dir if not specified
    if not pipeline_dir:
        pipeline_dir = env_defaults.get('pipeline-dir', pipeline_dir)

    # Add the environment-specific profile to the default profile
    if 'profile' in env_defaults:
        profile = f"{profile},{env_defaults['profile']}"

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

    # Print all parameters
    print("Running gOS with the following parameters:")
    print(f"Pipeline directory: {pipeline_dir}")
    print(f"Parameters file: {params_file}")
    print(f"Profile: {profile}")
    print(f"Resume: {resume}")

    # Load required modules
    load_modules_command = load_required_modules(env_defaults)

    # Initialize runner
    runner = NextflowRunner()

    # Build the command string without 'module load nextflow &&'
    command = (
        f"{load_modules_command} "
        f"{runner.cmd} run {pipeline_dir} "
        f"-params-file {params_file} "
        f"-profile {profile} "
        f"-with-report report_{runner.get_timestamp()}.html "
        f"-with-trace"
    )

    if resume:
        command += " -resume"

    # Execute the command using the updated runner
    runner.run(command)
