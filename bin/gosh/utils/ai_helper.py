import re
from typing import List

def extract_error_messages(log_content: str) -> List[str]:
    """
    Extract error messages from a Nextflow log file.
    
    Args:
        log_content (str): The content of the nextflow.log file
        
    Returns:
        List[str]: List of extracted error messages
    """
    # Regular expressions for matching
    timestamp_pattern = r'^[A-Z][a-z]{2}-\d{2}\s\d{2}:\d{2}:\d{2}\.\d{3}'
    error_pattern = r'ERROR'
    
    error_messages = []
    current_error = []
    in_error_block = False
    
    # Split the log content into lines
    lines = log_content.split('\n')
    
    for line in lines:
        # Check if line starts with timestamp
        is_timestamp_line = bool(re.match(timestamp_pattern, line))
        
        # If we're in an error block and find a new timestamp, 
        # we've reached the end of the current error
        if in_error_block and is_timestamp_line:
            if current_error:  # Only add if we have content
                error_messages.append('\n'.join(current_error))
                current_error = []
            in_error_block = False
            
        # Check if this is the start of a new error block
        if is_timestamp_line and error_pattern in line:
            in_error_block = True
            current_error = [line]
        # If we're in an error block, keep adding lines
        elif in_error_block:
            current_error.append(line)
            
    # Don't forget to add the last error block if we have one
    if current_error:
        error_messages.append('\n'.join(current_error))
        
    return error_messages

def read_log_file(log_path: str) -> str:
    """
    Read the contents of a log file.
    
    Args:
        log_path (str): Path to the nextflow.log file
        
    Returns:
        str: Contents of the log file
    """
    try:
        with open(log_path, 'r') as f:
            return f.read()
    except Exception as e:
        raise Exception(f"Failed to read log file: {str(e)}")
