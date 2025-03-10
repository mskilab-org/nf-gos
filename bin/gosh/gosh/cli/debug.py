#!/usr/bin/env python3

from os.path import exists
from csv import DictWriter
from sys import stdout
import click

@click.group(name='debug')
def debug_cli():
    """Debug commands for analyzing Nextflow runs."""
    pass

@debug_cli.command()
@click.argument('log_file', type=click.Path(exists=True), required=False)
def eye(log_file):
    """Analyze Nextflow log file using AI assistance."""
    try:
        # If no log file specified, look for .nextflow.log in current directory
        if not log_file:
            default_log = '.nextflow.log'
            if exists(default_log):
                log_file = default_log
            else:
                click.secho("Error: No .nextflow.log file found in current directory.", fg='red')
                click.secho("Please specify a log file path or run from a directory containing .nextflow.log", fg='yellow')
                return

        # Read the log file
        with open(log_file, 'r') as f:
            log_content = f.read()

        from ..utils.ai_helper import extract_error_messages, get_error_analysis_and_solution
        # Extract error messages
        error_messages = extract_error_messages(log_content)

        if not error_messages:
            click.secho("No errors found in the log file.", fg='green')
            return

        # Get AI analysis and solution
        analysis = get_error_analysis_and_solution(error_messages)

        # Print the analysis with colors
        click.secho("\n=== AI Analysis of Nextflow Errors ===", fg='blue', bold=True)
        click.echo("\n" + analysis)

    except Exception as e:
        click.secho(f"Error analyzing log file: {str(e)}", fg='red')

@debug_cli.command()
@click.option('-s', '--sample_names', type=str, help='Comma-separated list of sample ID(s).')
@click.option('-p', '--process_names', type=str, help='Comma-separated list of process name(s).')
@click.option('-o', '--output', type=click.Path(), help='Output file to save the results as CSV.')
def log(sample_names, process_names, output):
    """Retrieve Nextflow log entries based on samples or processes."""

    from ..core.nextflow_log import (
        get_all_entries,
        get_entries_with_sample_names,
        get_entries_with_process_names
    )

    try:
        entries_sample = []
        entries_process = []

        if sample_names:
            sample_list = [s.strip() for s in sample_names.split(',')]
            click.secho(f"Retrieving entries for sample name(s): {', '.join(sample_list)}", fg='blue')
            entries_sample = get_entries_with_sample_names(sample_list)

        if process_names:
            process_list = [p.strip() for p in process_names.split(',')]
            click.secho(f"Retrieving entries for process name(s): {', '.join(process_list)}", fg='blue')
            entries_process = get_entries_with_process_names(process_list)

        if sample_names and process_names:
            # Convert lists of dicts to sets of tuples for efficient intersection
            set_entries_sample = set(tuple(entry.items()) for entry in entries_sample)
            set_entries_process = set(tuple(entry.items()) for entry in entries_process)

            # Compute the intersection
            intersection = set_entries_sample & set_entries_process

            # Convert back to list of dicts
            entries = [dict(entry) for entry in intersection]
        elif sample_names:
            entries = entries_sample
        elif process_names:
            entries = entries_process
        else:
            click.secho("Retrieving all entries.", fg='blue')
            entries = get_all_entries()

        if not entries:
            click.secho("No log entries found.", fg='yellow')
            return

        # Remove duplicate entries if both sample and process names overlap
        unique_entries = [dict(t) for t in {tuple(d.items()) for d in entries}]

        if output:
            with open(output, 'w', newline='') as csvfile:
                writer = DictWriter(csvfile, fieldnames=unique_entries[0].keys())
                writer.writeheader()
                writer.writerows(unique_entries)
            click.secho(f"Results saved to {output}", fg='green')
        else:
            writer = DictWriter(stdout, fieldnames=unique_entries[0].keys())
            writer.writeheader()
            writer.writerows(unique_entries)

    except Exception as e:
        click.secho(f"Error retrieving log entries: {str(e)}", fg='red')
