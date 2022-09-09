#! /home/omelse/anaconda3/envs/AlphaFoldScripts/bin/python3
#Author: Okke Melse
#Last modified: 2022-02-09
import random
import string
import sys
import os
import dash
import io, base64
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import Bio.PDB as biopdb
from dash import Dash, dcc, Output, Input, State, html
import dash_bootstrap_components as dbc

def process_pdbfiles(pdb_file, session_id):
    content_type, content_pdb = pdb_file.split(',')
    decoded_pdb = base64.b64decode(content_pdb)
    with open(f"tmpFiles/{session_id}-AFpdb.pdb", "w") as f:
        f.write(decoded_pdb.decode("utf-8"))
    return 0

def get_IDDT(fileName, session_id):
    # Get IDDT for each residue (saved in B-factor)
    residueNrs = []
    IDDT_values = []
    pdbparser = biopdb.PDBParser()
    process_pdbfiles(fileName, session_id)
    pdb_protein = pdbparser.get_structure("protein", f"tmpFiles/{session_id}-AFpdb.pdb")

    # Only analyze first chain, remove others if present
    if len(list(pdb_protein.get_chains())) > 1:
        sys.stderr.write("More chains found, only chain A analyzed.\n")
        chainIDs = [x.get_id() for x in list(pdb_protein.get_chains())]
        for chainID in chainIDs[1:]:
            list(pdb_protein.get_models())[0].detach_child(chainID)

    # Collect IDDT values (b-factor)
    for residue in pdb_protein.get_residues():
        residueNrs.append(residue.get_id()[1])
        IDDT_values.append(residue["CA"].get_bfactor())
    return residueNrs, IDDT_values

## DASHBOARD ##
app = Dash(external_stylesheets=([dbc.themes.COSMO]))
dash.register_page(__name__,
                   title="Alphafold IDDT plot",
                   name="Generate Alphafold IDDT plot")

layout = dbc.Container([

    html.Div(
        html.H1("Plot confidentiality score from AlphaFold structures", style={'ext-align': 'center'})),

    html.Br(),

    html.Div([
        html.P(["Enter here the PDB file with IDDT values in B-factor column. Processing will start automatically after upload, please not this may take a couple of seconds.", html.Br()]),
        dcc.Upload(id="upload-pdb",
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={'width' : '100%',
                   'height' : '60px',
                   'lineHeight' : '60px',
                   'borderWidth': '1px',
                   'borderStyle': 'dashed',
                   'borderRadius': '5px',
                   'textAlign': 'center',
                   'margin': '10px'
                },
            multiple=False
        ),
    ]),

    # Show plot
    html.Div([
    html.Br(),
    html.Img(id="show_plot"),
    html.Br(),
    ])
    ])

@dash.callback(
    [Output(component_id="show_plot", component_property="src")],
    # Output(component_id="output_container2", component_property="children")],
    #[Input(component_id="submit-button-AF", component_property="n_clicks")],
    #[State(component_id="pdb_input", component_property="value")]
    [Input(component_id='upload-pdb', component_property="contents")],
    prevent_initial_call=True,
)
def make_plot(pdb_inputs):
    session_id = (''.join(random.choice(string.ascii_lowercase) for i in range(10)))

    ## matplotlib settings
    axistitleFont = {'fontname':'Arial', 'fontsize':12, 'fontweight':'normal'} #x-labels (set..)
    plt.rc('ytick', labelsize=12)
    plt.rc('xtick', labelsize=12)

    #Plotting
    fig = plt.figure(facecolor="white", edgecolor="white", figsize=(7,3))
    ax = fig.add_subplot(111)
    buf = io.BytesIO()

    p = []
    resNRS, IDDT = get_IDDT(pdb_inputs, session_id)
    p.append(ax.plot(resNRS, IDDT, linewidth=1.0))
    ax.set_xlabel("ResidueNr", **axistitleFont)
    ax.set_ylabel("IDDT", **axistitleFont)
    ax.set_ylim(0,100)
    plt.tight_layout()
    plt.savefig(buf, format="png")
    plt.close()
    data = base64.b64encode(buf.getbuffer()).decode("utf8")
    for file in os.scandir('tmpFiles'):
        if file.name.startswith(session_id):
            os.remove(file)
    return ["data:image/png;base64,{}".format(data)]
