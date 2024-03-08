def trim_content(content, x, y):
    if isinstance(content, str):
        return content[x : -y or None]
    elif isinstance(content, list):
        return content[x : -y or len(content)]
    else:
        return "Invalid content type. Please provide a string or a list."


def test_trim_content():
    seq = "ACGT"
    print(trim_content(seq, 0, 0))
