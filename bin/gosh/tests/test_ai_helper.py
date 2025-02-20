#!/usr/bin/env python3

import os
import sys
import click

# Add the bin directory to the Python path so we can import the module
sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from gosh.utils.ai_helper import extract_error_messages, query_ai, get_error_analysis_and_solution

def test_error_extraction():
    """Test the error message extraction functionality"""
    # Sample log content with an error message
    sample_log = """Dec-03 12:39:54.504 [Task monitor] INFO Some regular message
Dec-03 12:39:54.556 [Task monitor] ERROR nextflow.processor.TaskProcessor - Error executing process
Some error details
More error details
Dec-03 12:39:55.000 [Task monitor] INFO Another regular message"""

    # Extract error messages
    errors_message = extract_error_messages(sample_log)

    # Basic assertions
    assert errors_message, "Expected error messages"

    click.echo("✓ Error extraction test passed")
    return True

def test_error_extraction_from_logs():
    """Test error extraction from various Nextflow log files"""

    # Get path to test data directory
    test_dir = os.path.join(os.path.dirname(__file__), "test_data")

    # Test log file with no errors
    working_log = os.path.join(test_dir, "nextflow_working_heme_run.log")
    with open(working_log) as f:
        log_content = f.read()
    errors = extract_error_messages(log_content)
    assert len(errors) == 0, "Should find no errors in working log"

    # Test error logs
    error_logs = [
        "nextflow_failed_heme_run.log",
        "nextflow_failed_small_run.log",
        "nextflow_failed_full_run.log"
    ]

    for log_file in error_logs:
        log_path = os.path.join(test_dir, log_file)
        with open(log_path) as f:
            log_content = f.read()
        errors = extract_error_messages(log_content)
        print(log_path, '\n')
        print(errors, '\n\n')
        assert len(errors) > 0, f"Should find errors in {log_file}"

def test_query_ai():
    """Test querying the AI for error messages"""
    response = query_ai("tell me a joke")
    print(response)
    assert response, "Expected a response from the AI"

    click.echo("✓ Basic AI query test passed")


def test_error_analysis():
    """Test the error analysis and solution generation"""
    test_dir = os.path.join(os.path.dirname(__file__), "test_data")
    # use the working log file (should return a ValueError because of empty string)
    working_log = os.path.join(test_dir, "nextflow_working_heme_run.log")
    with open(working_log) as f:
        log_content = f.read()
    errors = extract_error_messages(log_content)

    # use the error log files
    error_logs = [
        "nextflow_failed_heme_run.log",
        "nextflow_failed_small_run.log",
        "nextflow_failed_full_run.log"
    ]

    for log_file in error_logs:
        log_path = os.path.join(test_dir, log_file)
        with open(log_path) as f:
            log_content = f.read()
        errors = extract_error_messages(log_content)
        analysis = get_error_analysis_and_solution(errors)
        print(log_file, '\n')
        print(analysis, '\n\n')
        assert analysis, f"Expected analysis for {log_file}"

    click.echo("✓ Error analysis test passed")


@click.command()
def main():
    """Run all tests for ai_helper.py"""
    click.echo("Running ai_helper.py tests...")

    try:
        # test_error_extraction()
        # test_error_extraction_from_logs()
        # test_query_ai()
        test_error_analysis()
        click.echo("\nAll tests passed! 🎉")
    except AssertionError as e:
        click.echo(f"\n❌ Test failed: {str(e)}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"\n❌ Unexpected error: {str(e)}", err=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
