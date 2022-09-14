#! /home/omelse/anaconda3/envs/EnvPy38/bin/python3
#Author: Okke Melse
#Last modified: 2022-09-08
import sys
import random
import string
from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import dash
from dash import Dash, dcc, Output, Input, State, html, ctx, dash_table
import dash_bootstrap_components as dbc
from utils.fileutils import save_multiple_files, empty_tmpFiles
import pandas as pd

# This scripts trims sequences and creates consensus sequences (fwd, fwd; or fwd, rev)
def make_consensus_seq(algn, seq1: SeqRecord, seq2: SeqRecord):
    """
    Takes pairwise2 alignment as input, returns consensus sequence
    :param algn: pairwise2 alignment
    :param seq1: SecRecord of first sequence (to retrieve phred quality data from)
    :param seq2: SecRecord of second sequence (to retrieve phred quality data from)
    :type algn: pairwise2.Alignment
    :type seq1: SeqRecord
    :type seq2: SeqRecord
    :return consensus sequence
    """
    seq1nr: int = 0 # to keep track of the nucl. sequence of this sequence in order to get phred quality
    seq2nr: int = 0 # to keep track of the nucl. sequence of this sequence in order to get phred quality
    consensus = []
    for nuc1, nuc2 in zip(algn[0], algn[1]):
        if nuc1 == nuc2:
            consensus.append(nuc1)
            seq1nr += 1
            seq2nr += 1
        elif nuc1 == '-':
            consensus.append(nuc2)
            seq2nr += 1
        elif nuc2 == '-':
            consensus.append(nuc1)
            seq1nr += 1
        elif nuc1 == 'N':
            consensus.append(nuc2)
            seq1nr += 1
            seq2nr += 1
        elif nuc2 == 'N':
            consensus.append(nuc1)
            seq1nr += 1
            seq2nr += 1
        elif seq1.letter_annotations['phred_quality'][seq1nr] > seq2.letter_annotations['phred_quality'][seq2nr]: #mismatch, take the base with highest phred quality
            consensus.append(nuc1)
            seq1nr += 1
            seq2nr += 1
        else:
            consensus.append(nuc2)
            seq1nr += 1
            seq2nr += 1

    consensus_record: SeqRecord = SeqRecord(Seq(''.join(consensus)),
                                                id=f"Consensus_{seq1.id}_{seq2.id}",
                                                description=f"{seq1.id} and {seq2.id}")
    return consensus_record

def parse_ab1_v2(filename1, filename2, direction="ff"):
    """
    Parses AB1 files with AbiIterator, thereby conserving phred quality data
    :param filename1: filename of sequence 1
    :param filename2: filename of sequence 2
    :param direction: string, either 'ff' or 'fr', where in case of the second, the second sequence (filename2) is reverse complemented
    """
    seq_adi_1: SeqIO.AbiIO.AbiIterator = SeqIO.AbiIO.AbiIterator(filename1, trim=True)
    seq_adi_2: SeqIO.AbiIO.AbiIterator = SeqIO.AbiIO.AbiIterator(filename2, trim=True)
    seq1: SeqRecord = next(seq_adi_1)
    seq1.id = ''.join(filename1.replace(" ","_").split("-")[1:])
    seq2: SeqRecord = next(seq_adi_2)
    seq2.id = ''.join(filename2.replace(" ","_").split("-")[1:])
    if direction == "fr":
        seq2 = seq2.reverse_complement(id=f"reverse_complement_of_{seq2.id}")
    return seq1, seq2

def get_basename(filenames):
    """
    Gets the 'basename' of each file (i.e. everything before '-') and saves in dictionary
    :param filenames: list of filenames to parse
    returns {"basename1": filename1, "basename2": filename2, ...}
    """
    file_basenames = {}
    for filename in filenames:
        # Check if files have expected .ab1 extension
        if not filename.split('.')[-1] == "ab1":
            raise RuntimeError(f"ERROR: file '{filename}' does not have the expected .ab1 extension")
        else:
            # Add basename and file in dict if not already present
            key = filename.split('-')[0]
            if key in file_basenames.keys():
                raise RuntimeError(f"ERROR: found multiple occurrences of file with basename {key}")
            else:
                file_basenames[filename.split('-')[0]] = filename
    return file_basenames

## DASHBOARD ##

app = Dash(external_stylesheets=([dbc.themes.ZEPHYR]))
dash.register_page(__name__,
                   title="Create consensus sequences",
                   name="Trimming and consensus sequences")

layout = dbc.Container([

    html.Div(
        html.H1("Calculate consensus sequences", style={'ext-align': 'center'})),

    html.Br(),

    html.Div([
        html.P([
            "This scipt will combine forward and reverse (or twice forward) sequencing results, "
            "trims the low-quality regions, and calculates the consensus sequence.\n"
            "The latter is performed as follows:\n"
            "1. Reverse complement is created if necessary\n"
            "2. Two sequences are aligned (local), with high gap-penalty\n"
            "3. Consensus sequence is calculated:\n"
            "   a. In non-aligned regions, the coding parent is selected\n"
            "   b. If one sequence contains a 'N', the other sequence is selected\n"
            "   c. In case of a mismatch, the sequence (nucleotide) with highest phred quality is selected\n"
        ],
        style={"whiteSpace": "pre"}),
        html.P([
            "IMPORTANT: Use the '-' character in your filenames to discriminate forward1/forward2, or forward/reverse.",
        ],
        style={"color": "red"}),
        html.P([
            "Thus the part before '-' has to be identical between sequencing data for the same gene, but unique"
            " (different) for sequencing data of other genes.\n"
            "For example: 'chimera1111-fw.ab1' and 'chimera1111-rv.ab1' will be recognized to belong together.\n"
        ],
        style={"whiteSpace": "pre"}),
        html.Br()
    ]),

    html.Div([
        html.P("Select direction of your sequences:"),
        dcc.Dropdown(id='direction',
                     options=[{"label": "Forward-Forward", "value": "ff"},
                              {"label": "Forward-Reverse", "value": "fr"}],
                     value="fr",
                     className="d-grid gap-2",
                     clearable=False,
                     style={"width": "50%"}
                     ),
        html.Br(),
    ]),

    html.Div([
        dcc.Upload(id="upload_seqs1",
                   children=html.Div([
                       "DNA forward sequences 1 in ABI format    (Drag and Drop or click to Select Files)"
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

        html.Em(id="out_filename1", className="text-info"),

        html.Br(),

        dcc.Upload(id="upload_seqs2",
                   children=html.Div([
                       "DNA reverse sequences (or forward sequences 2) in ABI format    (Drag and Drop or click to Select Files)"
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

        html.Em(id="out_filename2", className="text-info"),

    ]),

    html.Div([
        html.Br(),
        dbc.Button("Submit and calculate consensus",
                   color="primary",
                   id="submit-button-consensus",
                   n_clicks=0,
                   className="d-grid gap-2")
    ]),

    html.Div([dcc.Download(id="downloadoverview")]),
    html.Div([dcc.Download(id="downloadcconsensus")]),

    html.Br(),
    html.Div(id="output_msg_consensus",
             style={"whiteSpace": "pre"}),

])
# Callback to show uploaded files
@dash.callback(
    [Output(component_id="out_filename1", component_property="children"),
    Output(component_id="out_filename2", component_property="children")],
    [Input(component_id="upload_seqs1", component_property="filename"),
     Input(component_id="upload_seqs2", component_property="filename")],
    # prevent_initial_call=True,
)
def update_filenames(filenames1=None, filenames2=None):
    return [', '.join(filenames1) if filenames1 else "",
            ', '.join(filenames2) if filenames2 else ""]

# Callback to perform analysis
@dash.callback(
    [Output(component_id='downloadoverview', component_property='data'),
     Output(component_id='downloadcconsensus', component_property='data'),
     Output(component_id='output_msg_consensus', component_property='children'),
     Output(component_id='output_msg_consensus', component_property='className')],
    [Input(component_id='submit-button-consensus', component_property='n_clicks')],
    [State(component_id='direction', component_property='value'),
     State(component_id='upload_seqs1', component_property='filename'),
     State(component_id='upload_seqs1', component_property='contents'),
     State(component_id='upload_seqs2', component_property='filename'),
     State(component_id='upload_seqs2', component_property='contents')],
    prevent_initial_call=True,
)
def trim_and_consensus(n_clicks, direction, seq1filenames, seq1filedatas, seq2filenames, seq2filedatas):
    session_id = (''.join(random.choice(string.ascii_lowercase) for i in range(10)))

    # Step 1: Save files to disk
    try:
        save_multiple_files(seq1filedatas, seq1filenames, session_id, filetype='binary')
        save_multiple_files(seq2filedatas, seq2filenames, session_id, filetype='binary')
    except:
        empty_tmpFiles(session_id)
        return ["", "", "Some files could not be saved to disk, please check if you uploaded .ab1 files.", "text-danger"]

    # Step 2a: collect sequences which belong to each-other
    try:
        seq1_files: dict = get_basename(seq1filenames)
        seq2_files: dict = get_basename(seq2filenames)
    except RuntimeError as error:
        empty_tmpFiles(session_id)
        return ["", "", error.args[0], "text-danger"]

    # Step 2b: check if all files have a counterpart (i.e. reverse or second fwd seq)
    for key in seq1_files.keys():
        if not key in seq2_files.keys():
            empty_tmpFiles(session_id)
            return ["", "", f"ERROR: File with basename {key} not found in reverse/forward2 files.", "text-danger"]
    for key in seq2_files.keys():
        if not key in seq1_files.keys():
            empty_tmpFiles(session_id)
            return ["", "", f"ERROR: File with basename {key} not found in forward1 files.", "text-danger"]

    # Step 3: parse all ab1 files
    consensus: list[list[str, str, str, str, str]]
    consensus = []  # lists all complement sequences

    with open(f"consensus_seqs.fasta", "w") as f:
        for basename in seq1_files.keys():
            seq1, seq2 = parse_ab1_v2(f"tmpFiles/{session_id}-{seq1_files[basename]}",
                                      f"tmpFiles/{session_id}-{seq2_files[basename]}",
                                      direction)
            alignment = pairwise2.align.localms(seq1.seq, seq2.seq, 2, -1, -50, -5)
            warning = "Unreliable alignment" if len(alignment) > 1 else "No warning"
            consensus.append([seq1_files[basename],
                                seq2_files[basename],
                                warning,
                                str(make_consensus_seq(alignment[0], seq1, seq2).seq),
                                pairwise2.format_alignment(*alignment[0], full_sequences=True)])
            f.write(make_consensus_seq(alignment[0], seq1, seq2).format("fasta"))

    df = pd.DataFrame(consensus,
                      columns=["Seq1", "Seq2", "Warnings", "Consensus", "Alignment"])

    empty_tmpFiles(session_id)

    return [dict(content=df.to_csv(index=False), filename="consensus_overview.csv"),
            dcc.send_file("consensus_seqs.fasta"),
            "Consensus sequences successfully computed, download should start shortly.",
            "text-success"]
