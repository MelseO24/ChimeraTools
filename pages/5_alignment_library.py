#! /home/omelse/anaconda3/envs/EnvPy38/bin/python3
#Author: Okke Melse
#Last modified: 2022-09-06
import random
import string
import copy
import re
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import dash
from dash import Dash, dcc, Output, Input, State, html, ctx, dash_table
import dash_bootstrap_components as dbc
import matplotlib
import matplotlib.pyplot as plt
from utils.sequtils import sequence_info
from utils.fileutils import save_file, save_multiple_files, empty_tmpFiles
import pandas as pd
matplotlib.use("Agg")

# This script aligns all sequenced genes to the library

def get_nr_of_gaps(alignment: str, seq1: Seq, seq2: Seq):
    """
    Calculates number of gaps in pairwise sequence alignment
    """
    sequence_length = min(len(seq1), len(seq2)) # the maximum number of matches possible
    nr_matches = alignment.count("|")
    nr_mismatches = alignment.count(".")
    nr_gaps = sequence_length - nr_matches - nr_mismatches
    return nr_gaps



def parse_ab1(filename, session_id="tmp"):
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
        newtrace, ORFs = _parse_single_ab1(name)
        trace.append(newtrace)
        for frame in ORFs:
            tmpfasta.append(frame)
    # write all sequences concenated in one fasta file
    with open(f"tmpFiles/{session_id}-collectedFastaFiles.fasta", "w") as f:
        for seq in tmpfasta:
            f.write(seq.format("fasta"))
    libraryseqs_prot = sequence_info(f"tmpFiles/{session_id}-collectedFastaFiles.fasta")
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
        record_copy.id = record.id + f"-orf{rf_counter+1}"
        translated_seqs.append(copy.deepcopy(record_copy))
    return trace, translated_seqs

def parse_fasta(filename, session_id="tmp", translate=False):
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

    with open(f"tmpFiles/{session_id}-collectedFastaFiles.fasta", "w") as fwrite:
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
                        record_copy.id = record.id + f"-orf{rf_counter+1}"
                        fwrite.write(record_copy.format("fasta"))
                else:
                    fwrite.write(record.format("fasta"))
    libraryseqs = sequence_info(f"tmpFiles/{session_id}-collectedFastaFiles.fasta")
    return libraryseqs

def plot_trace(trace, session_id="tmp"):
    # fig = plt.figure(facecolor="white", edgecolor="white", figsize=(7,3))
    # ax = fig.add_subplot(111)
    plt.plot(trace["CHANNEL1"][0:500], color="blue")
    plt.plot(trace["CHANNEL2"][0:500], color="red")
    plt.plot(trace["CHANNEL3"][0:500], color="green")
    plt.plot(trace["CHANNEL4"][0:500], color="yellow")
    plt.savefig(f"tmpFiles/{session_id}-plot.png")

## DASHBOARD ##

# Build your components
app = Dash(external_stylesheets=([dbc.themes.ZEPHYR]))
dash.register_page(__name__,
                   title="Align sequenced library",
                   name="Analyze library sequencing")

# Customize your own Layout
layout = dbc.Container([

    html.Div(
        html.H1("Analyze library sequencing results", style={'ext-align': 'center'})),

    html.Br(),

    html.Div([

        html.Div("Note that this alignment is based on the translated AA sequence"
                 " (translation of sequencing files is performed automatically."),

        html.Br(),

        html.P([
            "File format of sequencing DNA files:"]),
        dcc.Dropdown(options=['fasta', 'AB1'],
                     value='fasta',
                     placeholder="Select the file format",
                     id="sequencing_format",
                     className="d-grid gap-2",
                     style={"width": '50%'},
                     clearable=False
                     ),

        html.P([
            html.Br(),
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
        html.Em(id="filename1", className="text-info"),

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
        html.Em(id="filename2", className="text-info"),

        html.Br(),
        html.Br(),

        html.Div("Select alignment method, global is generally recommended, \n"
                 "'local' should only be used when the sequenced DNA is too short (e.g. only forward used)",
                 style={"whiteSpace": "pre"}),

        dcc.Dropdown(options=['global', 'local'],
                     value='global',
                     placeholder="Select alignment method",
                     id="alignment_method",
                     className="d-grid gap-2",
                     style={"width": '50%'},
                     clearable=False
                     ),
    ]),

    # submit button
    html.Div([
        html.Br(),
        html.Em("Calculation may take a while, as long as the active tab is called 'Updating...', the script is running."),
        dbc.Button("Submit and start alignment (be patient, may take a while)",
                   color="primary",
                   id="submit-button2",
                   n_clicks=0,
                   className="d-grid gap-2")
    ]),

        html.Div([dcc.Download(id="downloadAlignment")]),

        html.Br(),

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
     Output(component_id='output_container_seqAlign', component_property="children"),
     Output(component_id='output_container_seqAlign', component_property="className")],
    [Input(component_id='submit-button2', component_property='n_clicks')],
    [State(component_id='upload-ab1', component_property='contents'),
     State(component_id='upload-library', component_property='contents'),
     State(component_id='upload-ab1', component_property='filename'),
     State(component_id='upload-library', component_property='filename'),
     State(component_id='sequencing_format', component_property='value'),
     State(component_id='alignment_method', component_property='value')],
    prevent_initial_call=True,
)
def process_seqalignment(n_clicks, ab1files, libraryfile, ab1filenames, libraryfilename, seq_fmt, alignment_method):
    session_id = (''.join(random.choice(string.ascii_lowercase) for i in range(10)))

    filetype = "text" if seq_fmt == "fasta" else "binary"

    # Step 0: Some checks
    if None in [ab1files, libraryfile]:
        return ["", "ERROR: upload both requested files.", "text-danger"]

    # Step 1a: Save sequencing (DNA) ABI/FASTA file to disk
    try:
        _localfilenames_sequencing = save_multiple_files(ab1files, ab1filenames, session_id, filetype)
    except:
        empty_tmpFiles(session_id)
        return ["","Sequencing file could not be saved to disk, check if you selected the correct file type.", "text-danger"]

    # Step 1b: Parse sequencing (DNA) ABI/FASTA file
    seq_sequencing: sequence_info
    if seq_fmt == "AB1":
        try:
            traces, seq_sequencing = parse_ab1(_localfilenames_sequencing, session_id)
        except:
            empty_tmpFiles(session_id)
            return ["","Sequencing file (AB1) could not be translated, probably no start codon was found,\n"
                       "or some fasta sequences have identical names (see terminal output for a more detailed error message).", "text-danger"]
    else: #seq_fmt == "fasta"
        try:
            seq_sequencing = parse_fasta(_localfilenames_sequencing, session_id, translate=True)
        except:
            empty_tmpFiles(session_id)
            return ["","Sequencing file (fasta) could not be translated, probably no start codon was found, \n"
                       "or some fasta sequences have identical names (see terminal output for a more detailed error message).", "text-danger"]

    # Step 2: Save and parse library FASTA files to disk
    _localfilename_library = save_file(libraryfile, libraryfilename, session_id)
    seq_library: SeqIO
    seq_library = SeqIO.parse(_localfilename_library, "fasta")

    dataLibrary = []
    aligner = Align.PairwiseAligner(mode=alignment_method,
                                    match_score=2,
                                    mismatch_score=-1,
                                    target_internal_open_gap_score=-50,
                                    target_internal_extend_gap_score=-5,
                                    query_internal_open_gap_score=-50,
                                    query_internal_extend_gap_score=-5,
                                    target_left_open_gap_score=-50,
                                    target_left_extend_gap_score=-5,
                                    query_left_open_gap_score=-50,
                                    query_left_extend_gap_score=-5,
                                    )
    record: SeqRecord
    for record in seq_library:
        found_exact_match = False
        # Find exact match in library
        for seqid in seq_sequencing.get_all_ids():
            if seq_sequencing.get_seq(seqid) == record.seq:
                matched_seq = [seqid, "max", "-", "-", seq_sequencing.get_seq(seqid)]
                found_exact_match = True
                continue
        if not found_exact_match:
            # find closest aligned library sequence
            best_alignment = {"score": 0,
                              "best_target": "",
                              "alignment": "-"}
            for seqid in seq_sequencing.get_all_ids():
                alignments = aligner.align(record.seq, seq_sequencing.get_seq(seqid))
                if alignments[0].score > best_alignment["score"]:
                    best_alignment = {"score": alignments[0].score,
                                      "best_target": seqid,
                                      "alignment": str(alignments[0])}

            matched_seq = [f"No exact match, best aligned with '{best_alignment['best_target']}'",
                           best_alignment['score'],
                           best_alignment['alignment'].count('.'),
                           get_nr_of_gaps(str(best_alignment['alignment']),
                                              record.seq,
                                              seq_sequencing.get_seq(best_alignment['best_target'])),
                           best_alignment['alignment']]

        dataLibrary.append([record.id,
                            matched_seq[0],  # seqid of matched sequencing result
                            matched_seq[1],  # Alignment score
                            matched_seq[2],  # Nr of mismatches in alignment
                            matched_seq[3],  # Nr of gaps in alignment
                            record.seq,      # expected sequence from library
                            matched_seq[4]]) # translated sequence from sequencing result

    df = pd.DataFrame(dataLibrary, columns=["Chimera",
                                            "Matched sequencing id",
                                            "Alignment score",
                                            "Mismatches",
                                            "Gaps",
                                            "Expected sequence (library)",
                                            "Translated sequencing result or best alignment"])

    empty_tmpFiles(session_id)
    return [dict(content=df.to_csv(index=False), filename="LibrarySequencingAnalysis.csv"),
            "Alignment successful, result will be downloaded automatically",
            "text-success"]
