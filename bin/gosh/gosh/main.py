import click
from .cli.run import run_cli
from .cli.debug import debug_cli

@click.group()
def cli():
    """gOSh - gOS sHell"""
    pass

# Register command groups
cli.add_command(run_cli, name='run')
cli.add_command(debug_cli, name='debug')

if __name__ == '__main__':
    cli()
