import dash
from dash import Dash, html, dcc
import dash_bootstrap_components as dbc


app = Dash(__name__, use_pages=True, external_stylesheets=([dbc.themes.SOLAR]))

app.layout = dbc.Container([
    html.H1("GLYCOtools", style={'ext-align': 'center'}),
    html.Div([
        html.Div(
            dcc.Link(
                f"{page['name']}", href=page["relative_path"], style={"font-size" : "18px"}
            )
        )
        for page in dash.page_registry.values()
    ]),
    html.Br(),
    dash.page_container,
    html.Footer(id = "footer",
        children = "\nCopyright CBR - Technical University of Munich.\nLast update: September 6th, 2022.",
        style={"whiteSpace" : "pre"}
        ),
    html.A("CBR Website", href='https://www.cbr.cs.tum.de', target='blank')
])

# Run app
if __name__=='__main__':
    app.run_server(host='0.0.0.0')

#if __name__=='__main__':
#    app.run_server(port=8052, debug=False)
