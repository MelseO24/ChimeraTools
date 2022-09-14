#! /home/omelse/anaconda3/envs/EnvPy38/bin/python3
#Author: Okke Melse
#Last modified: 2022-07-06

import sys
import os
import random
import string
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import dash
from dash import Dash, dcc, Output, Input, State, html
import dash_bootstrap_components as dbc
from utils.fileutils import empty_tmpFiles

# This script creates chimeric sequences out of parent sequences

#REQUIRED INPUT:
#parent_seqs: the fasta files of the parent sequences
#crossover_pts: the sequence nr of the 'previous' sequence to switch to new sequence.
#Format: [ [parent1Start, parent1End], [parent2Start, parent2End], .., [parentNStart, parentNend]]
#If start is start of sequence, use 0
#If end is end of sequence, use None
#e.g. [[0,255], [257, None]]

#parent_seqs = ["UGT71C1.fasta", "UGT71C2.fasta"]
#crossover_pts = [[0,255], [255, None]]

## Functions and classes ##

class sequence_info:
    def __init__(self, input_fasta):
        """Saves information from fasta file
         Provide input: fasta file, containing one or more sequences"""
        self.seqcounter = 0
        self.input_fasta = input_fasta
        self.fasta_seq = SeqIO.parse(open(input_fasta), 'fasta')

        self.seqdata = {}
        for fasta in self.fasta_seq:
            self.seqcounter += 1
            if fasta.id in self.seqdata.keys():
                raise KeyError(f"FATAL: multiple sequences with id: '{fasta.id}' found!.")
                sys.stderr.write(f"FATAL: multiple sequences with id: '{fasta.id}' found!.\n")
                usage()
                sys.exit(1)
            self.seqdata[fasta.id] = fasta.seq
        
        print(f"Analyzed {self.seqcounter} sequences.")

    def get_seq(self, seqid):
        return self.seqdata[seqid]

    def get_length(self, seqid):
        return len(self.seqdata[seqid])

    def get_all_ids(self):
        return self.seqdata.keys()

def usage():
    print("")
    print("USAGE:")
    print("Provide in script all parent sequences at 'parent_seqs'")
    print("Provide crossover points at crossover_pts:\n"
    "Format: [ [parent1Start, parent1End], [parent2Start, parent2End], .., [parentNStart, parentNend]]" \
    "If start is start of sequence, use '0'" \
    "If end is end of sequence, use ''" \
    "e.g. [[0,255], [257, ]]")
    return "Invalid"
        
def check_input(parent_seqs, crossover_pts):
    if len(parent_seqs) < 2:
        sys.stdout.write("FATAL: parent_seqs requires at least two parent fasta files\n")
        usage()
        sys.exit(1)
    elif len(crossover_pts) != len(parent_seqs):
        sys.stderr.write("FATAL: Number of provided parent sequences should equal number of crossover points provided.\n")
        usage()
        sys.exit(1)

    for filename in parent_seqs:
        if not os.path.isfile(filename):
            sys.stderr.write(f"FATAL, {filename} not found\n")
            usage()
            sys.exit(1)
    return 0

def combine_fasta_files(fasta_files, fasta_out):
    """ Reads in separate fasta files, and creates a single file containing all sequences"""
    sequences = []
    for filename in fasta_files:
        _seq = SeqIO.parse(open(filename), 'fasta')
        for fasta in _seq:
            sequences.append(fasta)
    SeqIO.write(sequences, fasta_out, 'fasta')
    return 0

def write_file(filename, value):
    f = open(filename, 'w')
    f.write(value)
    f.close()
    return 0

## Main ##
def create_chimeric_sequence(parent_seq1, parent_seq2, crossover_pts, session_id="tmp"):

    parent_seqs = [f"tmpFiles/{session_id}-seq1.fasta", f"tmpFiles/{session_id}-seq2.fasta"]

    write_file(parent_seqs[0], parent_seq1)
    write_file(parent_seqs[1], parent_seq2)
    
    # check_input(parent_seqs, crossover_pts)
    crossover_labels = [] #copy of crossover_pts with None replaced by 'end', for fasta header
    for cross_combi in crossover_pts:
        crossover_labels.append(["end" if value == None else value for value in cross_combi])

    # chimericSeq = The final sequence of the chimeric protein
    chimericSeq = SeqRecord(Seq(""),
                            id="chimeric_protein",
                            description=f"Parents: {' '.join(str(e) for e in parent_seqs)} -- crossoverPts: {','.join(str(e) for e in crossover_labels)}"
                            )

    # Create temporary file containing all parent sequences
    combine_fasta_files(parent_seqs, f"tmpFiles/{session_id}-tmp-fasta.fasta")

    # Create chimeric protein
    seqInfo = sequence_info(f"tmpFiles/{session_id}-tmp-fasta.fasta")
    for seqid, crossoverPt in zip(seqInfo.get_all_ids(), crossover_pts):
        chimericSeq.seq += seqInfo.get_seq(seqid)[crossoverPt[0] : crossoverPt[1]]
    
    ## Write final chimeric sequence to file
    #SeqIO.write(chimericSeq, 'chimericSeq.fasta', 'fasta')
    
    # print("Done: chimeric sequence saved in file: chimericSeq.fasta")
    return str(chimericSeq.format("fasta"))



## DASHBOARD ##

# Build your components
app = Dash(external_stylesheets=([dbc.themes.ZEPHYR]))
dash.register_page(__name__,
                   title="crossover_single",
                   name="In silico crossover - single")

# Customize your own Layout
layout = dbc.Container(
    html.Div([
        html.H1("Create chimeric sequences", style={'ext-align': 'center'}),
    
    html.Div([
        html.P(["Enter fasta 1 here: ", html.Br()]),
        dcc.Textarea(id='fasta1',
                  placeholder="Enter here your fasta of sequence 1",
                  style={'width' : '80%', 'height' : 150}),
        html.P(["Crossover position sequence 1: ", html.Br()]),
        dcc.Input(id="cut1",
                  placeholder='residuenr',
                  type="number",
                  size="20",
                  required=True
                  )
    ]),
    
    html.Br(),
    
   
    html.Div([
        html.P(["Enter fasta 2 here: ", html.Br()]),
        dcc.Textarea(id='fasta2',
                  placeholder="Enter here your fasta of sequence 2",
                  style={'width' : '80%', 'height' : 150}),
        html.P(["Crossover position sequence 2: ", html.Br()]),
        dcc.Input(id="cut2",
                  placeholder='residuenr',
                  type="number",
                  size="20",
                  required=True
                  )
    ]),
    
    html.Br(),
    
    dbc.Button("Submit and download chimeric fasta file",
               color="primary",
               id="submit-button",
               n_clicks=0,
               className="d-grid gap-2"),
    
    html.Br(),

    html.Div([dcc.Download(id="downloadFasta")]),
    
    html.Br(),
    
    html.Div(id='usedCuttingPoints'),
    
    #html.Br(),
    
    html.Div(id='output_container',
             style={"color" : "red"})

]))

# Callback allows components to interact
@dash.callback(
    # Output(component_id='output_container', component_property='children'),
    [Output(component_id="downloadFasta", component_property="data"),
     Output(component_id="output_container", component_property="children"),
     Output(component_id="usedCuttingPoints", component_property="children")],
    Input(component_id='submit-button', component_property='n_clicks'),
    [State(component_id='fasta1', component_property='value'),
     State(component_id='fasta2', component_property='value'),
     State(component_id='cut1', component_property='value'),
     State(component_id='cut2', component_property='value')],
    prevent_initial_call=True,
)
def calc_chimera(n_clicks, fasta1, fasta2, cut1, cut2):  # function arguments come from the component property of the Input
    session_id = (''.join(random.choice(string.ascii_lowercase) for i in range(10)))

    if n_clicks < 1:
        return ""
    elif (cut1 == None) or (cut2 == None):
        return ["", "Invalid, both crossover position 1 and 2 need to be provided", ""]
    elif not fasta1:
        return ["", "INVALID input, add a fasta sequence in Sequence 1", ""]
    elif not fasta1.startswith(">"):
        return ["", "INVALID input, sequence 1 is missing the fasta header '(>sequence)'", ""]
    elif not fasta2:
        return ["", "INVALID input, add a fasta sequence in Sequence 2", ""]
    elif not fasta2.startswith(">"):
        return ["", "INVALID input, sequence 2 is missing the fasta header '(>sequence)'", ""]
    else:
        crossover_pts = [[0,cut1], [cut2, None]]
        try:
            chimeric_seq = create_chimeric_sequence(fasta1, fasta2, crossover_pts, session_id)
        except KeyError as error:
            empty_tmpFiles(session_id)
            return ["", error.args[0], ""]

        empty_tmpFiles(session_id)

        return [dict(content=chimeric_seq,
                    filename="chimericprotein.fasta"),
                "Shuffling succeeded and download started",
                f"Used cutting points, fasta1: 0-{cut1}, fasta2: {cut2}-end."]

# Run app
if __name__=='__main__':
    app.run_server(port=8052, debug=False)


