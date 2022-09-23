import os
import base64

def _write_file(filename, value):
    with open(filename, "w") as f:
        f.write(value)

def save_file(file, filename, session_id="tmp", filetype='text'):
    """
    Saves uploaded file to disk
    :param file:
    :param filetype: 'text' or 'binary' (NOTE: ABI is binary!)
    :return: filename of fasta file on disk
    """
    content_type, content = file.split(',')
    decoded_file = base64.b64decode(content)
    if filetype == "text":
        with open(f"tmpFiles/{session_id}-{filename}", "w") as f:
            f.write(decoded_file.decode("utf-8"))
    elif filetype == "binary":
        with open(f"tmpFiles/{session_id}-{filename}", "wb") as f:
            f.write(decoded_file)
    return f"tmpFiles/{session_id}-{filename}"

def save_multiple_files(files, filenames, session_id="tmp", filetype='text'):
    """
    Saves uploaded files to disk
    :param filenames: list of files
    :param filetype: 'text' or 'binary' (NOTE: ABI is binary!)
    :return localfilenames: list of filenames on local disk
    """
    localfilenames = []
    for file, filename in zip(files, filenames):
        localfilenames.append(save_file(file, filename, session_id, filetype))
    return localfilenames

def empty_tmpFiles(session_id):
    for file in os.scandir('tmpFiles'):
        if file.name.startswith(session_id):
            try:
                os.remove(file)
            except:
                print(f"{file} could not be removed")