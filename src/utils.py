def get_file_lines(input_file_path: str) -> int:
    """
    Get the number of lines in a file.
    """
    with open(input_file_path, 'r') as f:
        return len(f.readlines())