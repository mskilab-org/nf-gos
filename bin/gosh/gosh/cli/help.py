from json import dump
import click

@click.group()
def help_cli():
    """Help command group for gOSh CLI."""
    pass

@help_cli.command()
@click.argument('query', type=str)
def ask(query):
    """Ask gosh a question about the nf-gOS pipeline"""
    from gosh.utils.ai_helper import answer_help_question, extract_new_params
    try:
        response = answer_help_question(query)
        click.echo("𓅃: " + response)
        # Try to extract new_params.json from the response
        try:
            new_params = extract_new_params(response)
            # Ask the user if they want to overwrite existing params.json
            if click.confirm("A new 'params.json' was provided in the response. Do you want to overwrite your existing 'params.json' with the new one?"):
                with open('params.json', 'w') as f:
                    dump(new_params, f, indent=4)
                click.echo("params.json has been updated.")
        except ValueError:
            # No new params.json provided; skip
            pass
    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
