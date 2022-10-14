# -*- coding: utf-8 -*-
#IMPORTANT: IN ORDER FOR THIS TO WORK, PYTHON2.7  NEEDS TO BE INSTALLED, AND THE SCHEME SCRIPTS NEED TO BE IN schema-tools (and executable)!
"""
Created on Aug 14

@author: Okke
"""

import random
import string
import dash
import subprocess
from dash import Dash, dcc, Output, Input, State, html
import dash_bootstrap_components as dbc
from utils.fileutils import save_file, empty_tmpFiles


## DASHBOARD ##

# Build components
app = Dash(external_stylesheets=([dbc.themes.ZEPHYR]))
dash.register_page(__name__,
                   title="SCHEMA",
                   name="SCHEMA analysis")

# Customize your own Layout
layout = dbc.Container(

    html.Div([
        html.Div([
            html.H1("Perform SCHEMA analysis", style={'ext-align': 'center'}),
            html.H5("See Saab-Rincon, Li, Meyer et al. Protein Engineering by Structure-Guided SCHEMA Recombination:"
                    " https://doi.org/10.1002/9783527634026.ch20", style={'ext-align': 'center'})
        ]),

    html.Br(),

    # PDB input
        html.P(["PDB file of parent 1 (either X-Ray or AlphaFold)",
                html.Br()]),

        dcc.Upload(id="input_pdb",
                   children=html.Div([
                       'Upload here the PDB file of parent 1 (Drag and Drop or click to Select Files)'
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
        html.Em(id="filename_input_pdb", className="text-info"),

        html.Br(),

        # MSA parents
        html.P(["Multiple sequence alignment (MSA) of all parents",
                html.Br()]),

        dcc.Upload(id="input_msa",
                   children=html.Div([
                       'Upload here the MSA in CLUSTAL format  (Drag and Drop or click to Select Files)'
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
        html.Em(id="filename_input_msa_parents", className="text-info"),

        html.Br(),

        # Alignment parent1 - PDB
        html.P(["Alignment between parent 1 and PDB sequence",
                html.Br()]),

        dcc.Upload(id="input_aln_parent_pdb",
                   children=html.Div([
                       'Upload here the alignment in CLUSTAL format  (Drag and Drop or click to Select Files)'
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
        html.Em(id="filename_input_aln_parent_PDB", className="text-info"),

        html.Br(),

        html.P([html.Br(), "Provide here the crossover points using the positions in the alignment (comma-separated):", html.Br()]),
        dcc.Input(id="input_crossPts",
                  placeholder="e.g. 68,209,259",
                  type="text",
                  size="30"
            ),
        html.Br(),
        html.Br(),

        # submit button
        html.Div([
            html.Br(),
            dbc.Button("Perform SCHEMA analysis",
                   color="primary",
                   id="input_submit",
                   n_clicks=0,
                   className="d-grid gap-2")
        ]),
        #Download output
        html.Div([dcc.Download(id="output_download")]),
        html.Br(),
        #Write comments
        html.Div(id='output_text',
                 style={"whiteSpace" : "pre"})

    ]),

)

@dash.callback(
    [Output(component_id='filename_input_pdb', component_property="children"),
     Output(component_id='filename_input_msa_parents', component_property="children"),
     Output(component_id='filename_input_aln_parent_PDB', component_property="children")],
    [Input(component_id='input_pdb', component_property='filename'),
     Input(component_id='input_msa', component_property='filename'),
     Input(component_id='input_aln_parent_pdb', component_property='filename')],
    prevent_initial_call=True,
)
def update_filenames(filename_input_pdb=None, filename_input_msa_parents=None, filename_input_aln_parent_PDB=None):
    # If upload triggered callback, show which file was uploaded, but do not start alignment.
    return [filename_input_pdb, filename_input_msa_parents, filename_input_aln_parent_PDB]

@dash.callback(
    [Output(component_id='output_download', component_property="data"),
     Output(component_id='output_text', component_property="children"),
     Output(component_id='output_text', component_property="className")],
    [Input(component_id='input_submit', component_property="n_clicks")],
    [State(component_id='input_pdb', component_property="contents"),
     State(component_id='input_msa', component_property="contents"),
     State(component_id='input_aln_parent_pdb', component_property="contents"),
     State(component_id='input_crossPts', component_property="value")],
    prevent_initial_call=True,
)
def run_schema(n_clicks, input_pdb, input_msa, input_aln_parent_pdb, input_crossPts):
    session_id = (''.join(random.choice(string.ascii_lowercase) for i in range(10)))

    if not input_crossPts:
        return ["", "ERROR, please provide crossover positions", "text-danger"]

    # save downloaded files to disk
    filenames = []
    for filename, data in zip(["input_pdb.pdb", "input_msa.aln", "input_aln_parent_pdb.aln"],
                              [input_pdb, input_msa, input_aln_parent_pdb]):
        filenames.append(save_file(data, filename, session_id, filetype="text"))

    # process crosspts
    input_crossPts = input_crossPts.split(",")
    input_crossPts = [int(x) for x in input_crossPts]
    input_crossPts.sort()
    input_crossPts = [str(x) for x in input_crossPts]
    with open(f"tmpFiles/{session_id}-xo.txt", "w") as f:
        f.write(" ".join(input_crossPts))

    # Run SCHEMA
    try:
        command_run_schemacontacts = (f"schema-tools/schemacontacts.py "
                                       f"-pdb {filenames[0]} "
                                       f"-msa {filenames[1]} "
                                       f"-pdbal {filenames[2]} "
                                       f"-o tmpFiles/{session_id}-contactmatrix.txt")

        process_schema_contacts = subprocess.Popen(command_run_schemacontacts,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True,
                                   universal_newlines=True)

        output_schema_contacts, _ = process_schema_contacts.communicate()

    except FileNotFoundError as error:
        return ["", f"ERROR, there is something wrong with your installation of Python2.7 or the SCHEMA scripts "
                    f"on your server.\nErrormsg: {error}",  "text-danger"]

    try:
        command_run_schemaenergy = (f"schema-tools/schemaenergy.py "
                                    f"-msa {filenames[1]} "
                                    f"-con tmpFiles/{session_id}-contactmatrix.txt "
                                    f"-xo tmpFiles/{session_id}-xo.txt "
                                     " -E"
                                     " -m")

        process_schemaenergy = subprocess.Popen(command_run_schemaenergy,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   shell=True,
                                   universal_newlines=True)

        output_schema_energy, _ = process_schemaenergy.communicate()

    except FileNotFoundError as error:
        return ["", f"ERROR, there is something wrong with your installation of Python2.7 or the SCHEMA scripts "
                    f"on your server.\nErrormsg: {error}",  "text-danger"]

    empty_tmpFiles(session_id)

    return [dict(content=output_schema_energy, filename="schema_output"),
            "Success",
            "text-success"]
