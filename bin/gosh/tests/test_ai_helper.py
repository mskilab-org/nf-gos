#!/usr/bin/env python3

import os
import sys
import click

# Add the bin directory to the Python path so we can import the module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from gosh.utils.ai_helper import extract_error_messages

def test_error_extraction():
    """Test the error message extraction functionality"""
    # Sample log content with an error message
    sample_log = """Dec-03 12:39:54.504 [Task monitor] INFO Some regular message
Dec-03 12:39:54.556 [Task monitor] ERROR nextflow.processor.TaskProcessor - Error executing process
Some error details
More error details
Dec-03 12:39:55.000 [Task monitor] INFO Another regular message"""

    # Extract error messages
    errors = extract_error_messages(sample_log)

    # Basic assertions
    assert len(errors) == 1, f"Expected 1 error message, got {len(errors)}"
    assert "ERROR" in errors[0], f"Expected ERROR in message, got: {errors[0]}"
    assert "Some error details" in errors[0], "Expected error details in message"

    click.echo("✓ Error extraction test passed")
    return True

@click.command()
def main():
    """Run all tests for ai_helper.py"""
    click.echo("Running ai_helper.py tests...")

    try:
        test_error_extraction()
        click.echo("\nAll tests passed! 🎉")
    except AssertionError as e:
        click.echo(f"\n❌ Test failed: {str(e)}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"\n❌ Unexpected error: {str(e)}", err=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
