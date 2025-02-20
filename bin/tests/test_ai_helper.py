import os
from gosh.utils.ai_helper import extract_error_messages

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
        "nextflow_error_1.log",
        "nextflow_error_2.log", 
        "nextflow_error_3.log"
    ]
    
    for log_file in error_logs:
        log_path = os.path.join(test_dir, log_file)
        with open(log_path) as f:
            log_content = f.read()
        errors = extract_error_messages(log_content)
        assert len(errors) > 0, f"Should find errors in {log_file}"
