#!/usr/bin/env python3

import os
import click
from utils.ai_helper import extract_error_messages, get_error_analysis_and_solution

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
            if os.path.exists(default_log):
                log_file = default_log
            else:
                click.secho("Error: No .nextflow.log file found in current directory.", fg='red')
                click.secho("Please specify a log file path or run from a directory containing .nextflow.log", fg='yellow')
                return

        # Read the log file
        with open(log_file, 'r') as f:
            log_content = f.read()

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
