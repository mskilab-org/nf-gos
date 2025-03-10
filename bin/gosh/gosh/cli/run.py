import os
import click
import socket
import shutil
import subprocess
import sys
from ..core.nextflow import NextflowRunner
from ..core.config import Config
from ..core.params_wizard import create_params_file
from ..core.nextflow_log import get_entries_with_process_names, get_entries_with_sample_names
from ..core.module_loader import get_environment_defaults, load_required_modules

@click.group(name='run')
def run_cli():
    """Run pipeline commands"""
    pass

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
@click.option('-p', '--processes',
              help='Comma-separated list of process names to rerun')
@click.option('-s', '--samples',
              help='Comma-separated list of sample IDs to rerun')
def pipeline(pipeline_dir, params_file, profile, resume, processes, samples):

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

    # Initialize a set to collect work directories to delete
    workdirs_to_delete = {}

    # Check if processes or samples are specified
    if processes or samples:
        # Parse the comma-separated lists into Python lists
        processes_list = [p.strip() for p in processes.split(',')] if processes else []
        samples_list = [s.strip() for s in samples.split(',')] if samples else []

        workdirs_processes = {}
        if processes_list:
            print("Retrieving work directories for specified processes...")
            entries = get_entries_with_process_names(processes_list)
            for entry in entries:
                workdir = entry['workdir']
                name = entry['name']
                if workdir not in workdirs_processes:
                    workdirs_processes[workdir] = set()
                workdirs_processes[workdir].add(name)

        workdirs_samples = {}
        if samples_list:
            print("Retrieving work directories for specified samples...")
            entries = get_entries_with_sample_names(samples_list)
            for entry in entries:
                workdir = entry['workdir']
                name = entry['name']
                if workdir not in workdirs_samples:
                    workdirs_samples[workdir] = set()
                workdirs_samples[workdir].add(name)

        if processes_list and samples_list:
            workdirs_set = set(workdirs_processes.keys()) & set(workdirs_samples.keys())
        elif processes_list:
            workdirs_set = set(workdirs_processes.keys())
        elif samples_list:
            workdirs_set = set(workdirs_samples.keys())
        else:
            workdirs_set = set()

        for workdir in workdirs_set:
            names = set()
            if workdir in workdirs_processes:
                names.update(workdirs_processes[workdir])
            if workdir in workdirs_samples:
                names.update(workdirs_samples[workdir])
            workdirs_to_delete[workdir] = names

        if workdirs_to_delete:
            print("The following work directories will be deleted:")
            print("Names, Work Directory")
            for workdir, names in sorted(workdirs_to_delete.items()):
                names_str = ', '.join(sorted(names))
                print(f"{names_str}, {workdir}")

            confirm = input("Do you want to delete these directories? (yes/no): ").strip().lower()
            if confirm == 'yes':
                confirm_again = input("Are you sure? This will cause the pipeline to rerun from these steps. (yes/no): ").strip().lower()
                if confirm_again == 'yes':
                    for workdir in workdirs_to_delete:
                        shutil.rmtree(workdir)
                    print("Directories deleted.")
                else:
                    print("Run cancelled.")
                    sys.exit(0)
            else:
                print("Run cancelled.")
                sys.exit(0)
        else:
            print("No matching work directories found to delete.")
            sys.exit(0)

    load_modules_command = load_required_modules(env_defaults)

    runner = NextflowRunner()

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

    runner.run(command)
