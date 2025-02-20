import re
from typing import List
from openai import OpenAI

ERROR_INTERPRETER_PROMPT = """You are a Nextflow pipeline debugging expert. Your task is to analyze the provided error messages and:
1. Identify the root cause of the error
2. Determine which part of the pipeline failed and for which sample name
3. Extract relevant details like file names/paths, samples, process names, commands, or parameters involved
4. Summarize this in a clear, technical manner
Be concise and focus only on the technical details of what went wrong.

Your response should be in the following format:

Sample Name:
Failed Process:
Work Directory:
Error:
Possible Causes:

"""

ERROR_ADVISOR_PROMPT = """You are a bioinformatics workflow expert specializing in Nextflow pipelines. Based on the error analysis provided:
1. Suggest specific steps to resolve the issue
2. Provide practical solutions that the user can implement
3. If relevant, explain any pipeline-specific considerations
4. If needed, recommend configuration changes or system requirements
Keep suggestions actionable and direct. Focus on practical solutions rather than theoretical explanations."""

def extract_error_messages(log_content: str) -> str:
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
    error_messages_str = ''
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

    for error in error_messages:
        error_messages_str += error + '\n\n'

    return error_messages_str

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

def query_ai(query: str, system_prompt: str = "You are a helpful assistant.") -> str:
    """
    Send a query to OpenAI API and get the response.

    Args:
        query (str): The user query to send to the AI
        system_prompt (str): The system prompt to set AI behavior/context

    Returns:
        str: The AI response text
    """
    if not query or query == "":
        raise ValueError("Query cannot be empty.")

    try:
        client = OpenAI()
        completion = client.chat.completions.create(
            model="gpt-4",
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": query}
            ]
        )
        return completion.choices[0].message.content
    except Exception as e:
        raise Exception(f"Failed to get AI response: {str(e)}")

def get_error_analysis_and_solution(error_messages: str) -> str:
    """
    Chain two AI queries to get both error interpretation and solution advice.

    Args:
        error_messages (str): The extracted error messages from the log

    Returns:
        str: Combined analysis and solution from the AI

    Raises:
        ValueError: If error_messages is empty
    """
    if not error_messages or error_messages.strip() == "":
        raise ValueError("No error messages found to analyze.")

    try:
        # First call: Interpret the error
        error_interpretation = query_ai(error_messages, ERROR_INTERPRETER_PROMPT)

        # Second call: Get solution based on interpretation
        solution = query_ai(error_interpretation, ERROR_ADVISOR_PROMPT)

        # Combine the responses with clear separation
        combined_response = f"""
        ERROR INTERPRETATION:\n{error_interpretation}\n\n
        RECOMMENDED SOLUTION:\n{solution}\n\n"""

        return combined_response

    except Exception as e:
        raise Exception(f"Failed to complete error analysis chain: {str(e)}")
