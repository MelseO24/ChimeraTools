#! /home/omelse/anaconda3/envs/EnvPy38/bin/python3
#Author: Okke Melse
#Last modified: 2022-09-06
import copy
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import base64
import dash
from dash import Dash, dcc, Output, Input, State, html, ctx, dash_table
import dash_bootstrap_components as dbc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from utils.sequtils import sequence_info
import pandas as pd

# This script aligns all sequenced genes to the library

def _write_file(filename, value):
    with open(filename, "w") as f:
        f.write(value)

def save_file(file, filename, filetype='text'):
    """
    Saves uploaded file to disk
    :param file:
    :param filetype: 'text' or 'binary' (NOTE: ABI is binary!)
    :return: filename of fasta file on disk
    """
    content_type, content = file.split(',')
    decoded_file = base64.b64decode(content)
    if filetype == "text":
        with open(f"tmp-{filename}", "w") as f:
            f.write(decoded_file.decode("utf-8"))
    elif filetype == "binary":
        with open(f"tmp-{filename}", "wb") as f:
            f.write(decoded_file)
    return f"tmp-{filename}"

def save_multiple_files(files, filenames, filetype='text'):
    """
    Saves uploaded files to disk
    :param filenames: list of files
    :param filetype: 'text' or 'binary' (NOTE: ABI is binary!)
    :return localfilenames: list of filenames on local disk
    """
    localfilenames = []
    for file, filename in zip(files, filenames):
        localfilenames.append(save_file(file, filename, filetype))
    return localfilenames

def parse_ab1(filename):
    """
    Parses ABI file(s) in biopython SeqIO record, and returns dictionary with channel readouts and FASTA
    :param filename: ABI filename (get from process_ab1)
    :return: list(dict(["CHANNEL1","CHANNEL2","CHANNEL3","CHANNEL4","FASTA_DNA","FASTA_AA"]))
    :return: sequence_info with parsed sequences
    """
    # Process trace file in biopython
    # trace = dict.fromkeys(["CHANNEL1","CHANNEL2","CHANNEL3","CHANNEL4"])
    trace = []
    tmpfasta = [] # collects all sequences in fasta format to later write in single file
    if type(filename) == str:
        filename = [filename]
    for name in filename:
        newtrace, readingFrames = _parse_single_ab1(name)
        trace.append(newtrace)
        for frame in readingFrames:
            tmpfasta.append(frame)
    # write all sequences concenated in one fasta file
    with open("tmp-collectedFastaFiles.fasta", "w") as f:
        for seq in tmpfasta:
            f.write(seq.format("fasta"))
    libraryseqs_prot = sequence_info("tmp-collectedFastaFiles.fasta")
    return trace, libraryseqs_prot

def _parse_single_ab1(filename):
    """
    Does the actual parsing of a ABI file, but users are recommended to use the wrapper parse_ab1 instead
    :param filename:
    :return: dict(["CHANNEL1","CHANNEL2","CHANNEL3","CHANNEL4","FASTA_DNA","FASTA_AA"])
    """
    trace = {}
    record: SeqRecord
    record = SeqIO.read(filename, "abi")
    for channelNr, data in enumerate(["DATA9","DATA10","DATA11","DATA12"]):
        trace[f"CHANNEL{channelNr+1}"] = record.annotations["abif_raw"][data]

    # Translate all reading frames
    translated_seqs = [] # list of seqrecords with all reading frames
    record_copy: SeqRecord
    record_copy = copy.deepcopy(record)
    for rf_counter, readingFrame in enumerate([m.start() for m in re.finditer("ATG", str(record.seq))]):
        record_copy.letter_annotations = {}
        record_copy.seq = record.seq[readingFrame:]
        record_copy.seq = record_copy.translate(to_stop=True).seq
        record_copy.id = record.id + f"-rf{rf_counter+1}"
        translated_seqs.append(copy.deepcopy(record_copy))
    return trace, translated_seqs

def parse_fasta(filename, translate=False):
    """
    Parses fasta file (with one or more sequences) in biopython SeqRecord
    :param filename: filename of fasta to parse
    :param translate: if input DNA needs to be translated to protein seq
    :type translate: bool
    :return: sequence_info with parsed sequences
    """
    libraryseqs: sequence_info
    # concentate file contents in a single file
    if type(filename) == str:
        filename = [filename]

    with open("tmp-collectedFastaFiles.fasta", "w") as fwrite:
        for name in filename:
            all_sequences = SeqIO.parse(name, "fasta")
            record: SeqRecord
            for record in all_sequences:
                if translate:
                    record_copy: SeqRecord
                    record_copy = copy.deepcopy(record)
                    for rf_counter, readingFrame in enumerate([m.start() for m in re.finditer("ATG", str(record.seq))]):
                        record_copy.seq = record.seq[readingFrame:]
                        record_copy.seq = record_copy.translate(to_stop=True).seq
                        record_copy.id = record.id + f"-rf{rf_counter+1}"
                        fwrite.write(record_copy.format("fasta"))
                else:
                    fwrite.write(record.format("fasta"))
    libraryseqs = sequence_info("tmp-collectedFastaFiles.fasta")
    return libraryseqs

def plot_trace(trace):
    # fig = plt.figure(facecolor="white", edgecolor="white", figsize=(7,3))
    # ax = fig.add_subplot(111)
    plt.plot(trace["CHANNEL1"][0:500], color="blue")
    plt.plot(trace["CHANNEL2"][0:500], color="red")
    plt.plot(trace["CHANNEL3"][0:500], color="green")
    plt.plot(trace["CHANNEL4"][0:500], color="yellow")
    plt.savefig("tmpplot.png")

## DASHBOARD ##

# Build your components
app = Dash(external_stylesheets=([dbc.themes.COSMO]))
dash.register_page(__name__,
                   title="Align sequenced library",
                   name="Analyze library sequencing")

# Customize your own Layout
layout = dbc.Container([
    html.Div([

        html.P([
                   "Upload here the sequencing DNA files (.fasta or .ab1), multiple files are allowed (keep ctrl pressed):\n",
                   html.Br()]),

        dcc.Upload(id="upload-ab1",
                   children=html.Div([
                       'DNA sequencing ABI/fasta File  (Drag and Drop or click to Select Files)'
                   ]),
                   style={'width': '100%',
                          'height': '60px',
                          'lineHeight': '60px',
                          'borderWidth': '1px',
                          'borderStyle': 'dashed',
                          'borderRadius': '5px',
                          'textAlign': 'center',
                          'margin': '10px'
                          },
                   multiple=True
                   ),
        html.Div(id="filename1",
                 style={"color":"white"}),
        html.Br(),
        dcc.Dropdown(options=['fasta', 'AB1'],
                     value='fasta',
                     placeholder="Select the file format",
                     id="sequencing_format",
                     className="d-grid gap-2",
                     style={"width": '80%'},
                     clearable=False
                     ),
        html.Br(),

        html.P([
            "Upload here the library protein sequences in FASTA format (single file):",
            html.Br()]),

        dcc.Upload(id="upload-library",
                   children=html.Div([
                       'Library amino-acid sequences  (Drag and Drop or click to Select Files)'
                   ]),
                   style={'width': '100%',
                          'height': '60px',
                          'lineHeight': '60px',
                          'borderWidth': '1px',
                          'borderStyle': 'dashed',
                          'borderRadius': '5px',
                          'textAlign': 'center',
                          'margin': '10px'
                          },
                   multiple=False
                   ),
        html.Div(id="filename2",
                 style={"color": "white"}),
        html.Br(),
    ]),

    # submit button
    html.Div([
        html.Br(),
        dbc.Button("Submit and start alignment",
                   color="primary",
                   id="submit-button2",
                   n_clicks=0,
                   className="d-grid gap-2")
    ]),

        html.Div([dcc.Download(id="downloadAlignment")]),

        html.Div(id="output_container_seqAlign",
                 style={"whiteSpace": "pre", "color": "red"}),

        html.Br(),

])

# Callback to show which file was uploaded
@dash.callback(
    [Output(component_id='filename1', component_property="children"),
     Output(component_id='filename2', component_property="children")],
    [Input(component_id='upload-ab1', component_property='filename'),
    Input(component_id='upload-library', component_property='filename')],
    prevent_initial_call=True,
)
def update_filenames(ab1filenames=None, libraryfilename=None):
    # If upload triggered callback, show which file was uploaded, but do not start alignment.
    return [', '.join(ab1filenames) if ab1filenames else "", libraryfilename]

# Callback to perform analysis
@dash.callback(
    [Output(component_id='downloadAlignment', component_property='data'),
     Output(component_id='output_container_seqAlign', component_property="children")],
    [Input(component_id='submit-button2', component_property='n_clicks')],
    [State(component_id='upload-ab1', component_property='contents'),
     State(component_id='upload-library', component_property='contents'),
     State(component_id='upload-ab1', component_property='filename'),
     State(component_id='upload-library', component_property='filename'),
     State(component_id='sequencing_format', component_property='value')],
    prevent_initial_call=True,
)
def process_seqalignment(n_clicks, ab1files, libraryfile, ab1filenames, libraryfilename, seq_fmt):

    filetype = "text" if seq_fmt == "fasta" else "binary"

    # Step 0: Some checks
    if None in [ab1files, libraryfile]:
        return ["ERROR: upload both requested files."]

    # Step 1a: Save sequencing (DNA) ABI/FASTA file to disk
    try:
        _localfilenames_sequencing = save_multiple_files(ab1files, ab1filenames, filetype)
    except:
        return ["","Sequencing file could not be saved to disk, check if you selected the correct file type."]

    # Step 1b: Parse sequencing (DNA) ABI/FASTA file to disk
    seq_sequencing: sequence_info
    if seq_fmt == "AB1":
        try:
            traces, seq_sequencing = parse_ab1(_localfilenames_sequencing)
        except:
            return ["","Sequencing file (AB1) could not be translated, probably no start codon was found."]
    else: #seq_fmt == "fasta"
        try:
            seq_sequencing = parse_fasta(_localfilenames_sequencing, translate=True)
        except:
            return ["","Sequencing file (fasta) could not be translated, probably no start codon was found."]

    # Step 2: Save and parse library FASTA files to disk
    _localfilename_library = save_file(libraryfile, libraryfilename)
    seq_library: SeqIO
    seq_library = SeqIO.parse(_localfilename_library, "fasta")

    dataLibrary = []
    record: SeqRecord
    for record in seq_library:
        matched_seq = ["No match", ""]  # [seqid, translated sequence]
        for seqid in seq_sequencing.get_all_ids():
            if seq_sequencing.get_seq(seqid) == record.seq:
                matched_seq = [seqid, seq_sequencing.get_seq(seqid)]
        dataLibrary.append([record.id,
                            matched_seq[0],  # seqid of matched sequencing result
                            record.seq,      # expected sequence from library
                            matched_seq[1]]) # translated sequence from sequencing result

    df = pd.DataFrame(dataLibrary, columns=["Chimera", "Matched sequencing id", "Expected sequence (library)", "Translated sequencing result"])

    return [dict(content=df.to_csv(index=False), filename="LibrarySequencingAnalysis.csv"),
            "Alignment successful, result will be downloaded automatically"]
