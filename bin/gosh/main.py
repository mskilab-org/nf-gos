import click
from cli.run import run_cli

@click.group()
def cli():
    """gOSh - gOS sHell"""
    pass

# Register command groups
cli.add_command(run_cli)

if __name__ == '__main__':
    cli()
