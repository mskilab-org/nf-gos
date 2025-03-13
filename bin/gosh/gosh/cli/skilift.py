from os import makedirs, path
from subprocess import run, CalledProcessError
import click

@click.group()
def skilift_cli():
    """Skilift command group for gOSh CLI."""
    pass

@skilift_cli.command()
@click.argument('results_dir', type=click.Path(), default='./results', required=False)
@click.option('-t', '--cohort_type', type=click.Choice(['paired', 'tumor_only', 'heme']), default='paired', help='Type of the cohort.')
@click.option('-o', '--gos_dir', type=click.Path(), required=True, help='Path to gos_dir.')
@click.option('-c', '--cores', type=int, default=1, help='Number of cores to use.')
def up(results_dir, cohort_type, gos_dir, cores):
    """Lift raw data into gOS compatible formats using Skilift."""
    if not path.isdir(results_dir):
        click.echo(f"Error: The directory '{results_dir}' does not exist.")
        return

    makedirs(expanduser(gos_dir), exist_ok=True)

    r_code = f'''
    devtools::load_all("~/git/skilift")
    cohort <- Cohort$new("{results_dir}", cohort_type="{cohort_type}")
    saveRDS(cohort, "{gos_dir}/cohort.rds")
    lift_all(cohort, output_data_dir="{gos_dir}", cores={cores})
    '''

    try:
        run(['Rscript', '-e', r_code], check=True)
    except CalledProcessError as e:
        click.echo(f"Error executing R script: {e}")
