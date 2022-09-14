import dash
from dash import Dash, html
import dash_bootstrap_components as dbc


app = Dash(__name__, use_pages=True, external_stylesheets=([dbc.themes.ZEPHYR]))

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "18rem",
    "padding": "2rem 1rem",
    # "background-color": "#f8f9fa",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "1rem 1rem"
}

sidebar = html.Div(
    [
        html.H2("GlycoTools", className="display-4"),
        html.Hr(),
        html.P(
            "Algorithms", className="lead"
        ),
        dbc.Nav(
            [
                dbc.NavLink(f"{page['name']}", href=page["relative_path"], style={"font-size": "18px"}) for page in dash.page_registry.values()
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)


# app.layout = html.Div([dcc.Location(id="url"), sidebar, content])
app.layout = dbc.Container([
    html.Br(),
    sidebar,
    html.Br(),
    html.Div([dash.page_container], style=CONTENT_STYLE),
    html.Footer(id = "footer",
        children = "\nCopyright CBR - Technical University of Munich.\nLast update: September 8th, 2022.",
        style={"whiteSpace" : "pre",
               "margin-left": "18rem",
               "margin-right": "2rem",
               "padding": "1rem 1rem"}
        ),
    html.A("CBR Website", href='https://www.cbr.cs.tum.de', target='blank', style=CONTENT_STYLE)
])

# Run app
if __name__=='__main__':
    app.run_server()

#if __name__=='__main__':
#    app.run_server(port=8052, debug=False)
