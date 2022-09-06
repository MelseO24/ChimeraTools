# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 15:22:10 2022

@author: Okke
"""

import sys
import os
import re
import itertools

import dash
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import dash
from dash import Dash, dcc, Output, Input, State, html
import dash_bootstrap_components as dbc

def write_file(filename, value):
    f = open(filename, 'w')
    f.write(value)
    f.close()
    return 0

def show_crossover_positions_in_alignment(msa, crossPts):
    """
    Returns a printable list of the sequences in MSA, focused on the residues
    found at the crossover positions

    Parameters
    ----------
    msa : SeqIO.SeqRecord
        Parsed MSA via AlignIO
    crossPts : list
        List of crossover positions (starts at 0).

    Returns
    -------
    Printable list of sequences in MSA at crossover positions.

    """
    printable_text = []
    for seq in msa:
        # printable_text += [str(seq.id), "\t"] + [str(seq.seq[x-1:x+2]) + "\t" for x in crossPts] + ["\n"]
        printable_text += [str(seq.id), "\t"] + [str(seq.seq[x-1:x+2]) + "\t" for x in crossPts] + ["\n"]
    return " ".join(printable_text)

def merge_sequences(sequences):
    """
    Merges sequences of list consisting of Bio.Seq sequences, and returns a string containing the concenated sequence
    Parameters
    ----------
    sequences: list of Bio.Seq

    Returns
    -------
    String of concenated sequences, with "-" removed
    """
    merged_seq = ""
    for fragment in sequences:
        merged_seq += str(fragment)
    merged_seq = re.sub('-', '', merged_seq)
    return merged_seq

def calc_library(msa, crossPts, label="full"):
    """
    Returns fasta files of all sequences resulting from shuffling

    Parameters
    ----------
    msa : SeqIO.SeqRecord
        Parsed MSA via AlignIO
    crossPts : list
        List of crossover positions (starts at 0).
        
    Returns
    -------
    Fasta file containing all sequences resulting from shuffling
    """

    # Phase 1: split sequences in fragments
    sequence_frags = []
    seqids = []
    for seq in msa:
        seqids.append(seq.id)
        sequence_frags.append([])
        for nr, cut in enumerate(crossPts):
            if nr == 0:
                sequence_frags[-1].append(seq.seq[:crossPts[nr]])
            elif nr == len(crossPts)-1:
                sequence_frags[-1].append(seq.seq[crossPts[nr-1]:crossPts[nr]])
                sequence_frags[-1].append(seq.seq[crossPts[nr]:])
            else:
                sequence_frags[-1].append(seq.seq[crossPts[nr-1]:crossPts[nr]])

    # Phase 2: create all possible chimeric sequences in library
    sequence_frags = np.array(sequence_frags, dtype=object).transpose() # Transpose to get [[frag1, frag1, frag1], [frag2, frag2, ]] etc.
    shuffled_fragments = itertools.product(*sequence_frags)
    shuffled_seq_index = itertools.product(*[list(range(len(seqids))) for x in range(len(crossPts)+1)]) #index values of shuffled domains, to retrieve from which id they came later on
    all_chimeric_seqs = []
    fasta_out = ""
    for sequence, seqindex in zip(shuffled_fragments, shuffled_seq_index):
        all_chimeric_seqs.append(SeqRecord(Seq(merge_sequences(sequence))))
        chimeraname = ''.join([str(x+1) for x in seqindex])
        if label == "full":
            label = ":".join([f"{chimeraname}:"
                              f"frag{fragNr+1}={seqids[seqind]}" for fragNr, seqind in enumerate(seqindex)])
        else:
            label = chimeraname
        all_chimeric_seqs[-1].id = label
        all_chimeric_seqs[-1].description = "Chimeric protein"
        fasta_out += all_chimeric_seqs[-1].format("fasta")
    return fasta_out

## DASHBOARD ##

# Build components
app = Dash(external_stylesheets=([dbc.themes.COSMO]))
dash.register_page(__name__,
                   title="Chimeric Library",
                   name="In silico library generation for chimeric proteins")

# Customize your own Layout
layout = dbc.Container(

    html.Div([
        html.Div(
            html.H1("Create libraries of chimeric proteins", style={'ext-align': 'center'})),
        # input MSA
        html.P(["Enter here the MSA (Clustal format):", html.Br()]),
        dcc.Textarea(id='msa_input',
                  placeholder="Enter here your MSA in a supported format",
                  style={'width' : '80%', 'height' : 250}),
        html.Br(),
        html.Br(),
        dcc.Dropdown(options=['Clustal','emboss','fasta'],
                     placeholder="Select an alignment format",
                     id="alignment_format",
                     className="d-grid gap-2",
                     style={"width" : '80%'}
                     ),
        # input crossoverPts
        html.P([html.Br(), "Provide here the crossover points using the positions in the alignment (comma-separated):", html.Br()]),
        dcc.Input(id="crossPts",
                  placeholder="e.g. 68,209,259",
                  type="text",
                  size="30"
            ),
        html.Br(),
        html.Br(),

        dcc.Dropdown(options=[
                    {"label": "Full annotation", "value": "full"},
                    {"label": "Short annotation", "value": "short"}],
                     id="annotation",
                     value="short",
                     className="d-grid gap-2",
                     clearable=False,
                     style={"width": '80%'}
                     ),

        # submit button
        html.Div([
            html.Br(),
            dbc.Button("Submit and download chimeric fasta file",
                   color="primary",
                   id="submit-button1",
                   n_clicks=0,
                   className="d-grid gap-2")
        ]),
        #Download fasta
        html.Div([dcc.Download(id="downloadFasta1")]),
        html.Br(),
        #Write cutting points
        html.Div(id='output_container1',
                 style={"whiteSpace" : "pre"})

    ]),

)
@dash.callback(
    [Output(component_id="downloadFasta1", component_property="data"),
     Output(component_id="output_container1", component_property="children")],
    [Input(component_id="submit-button1", component_property="n_clicks")],
    [State(component_id="msa_input", component_property="value"),
     State(component_id="crossPts", component_property="value"),
     State(component_id="alignment_format", component_property="value"),
     State(component_id="annotation", component_property="value")],
    prevent_initial_call=True,
)
def create_library(n_clicks, msa_input, crossPts_input, alignment_format, annotation):
    write_file("tmp-msa.aln", msa_input)
    try:
        align = AlignIO.read("tmp-msa.aln", alignment_format.lower())
    except ValueError:
        os.remove("tmp-msa.aln")
        return("", f"Error: MSA could not be parsed as {alignment_format} alignment. Are you sure you provided the correct format?\n")
    os.remove("tmp-msa.aln")
    crossPts_input = crossPts_input.split(",")
    crossPts = [int(x)-1 for x in crossPts_input]

    return [dict(content=calc_library(align, crossPts, annotation), filename="library_chimeras.fasta"),
            f"Calculation successful. Created {len(align) ** (len(crossPts)+1)} new sequences.\n"
            f"Below you find the shuffling positions (including one AA before and after the shuffling position).\n"
            f"Download of the library should start shortly.\n\n{show_crossover_positions_in_alignment(align, crossPts)}"]


# Run app
if __name__=='__main__':
    app.run(host='0.0.0.0')
    #app.run_server(port=8052, debug=False)
