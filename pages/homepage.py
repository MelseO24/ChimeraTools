import dash
from dash import Dash, dcc, Output, Input, State, html
import dash_bootstrap_components as dbc

# Build your components
app = Dash(external_stylesheets=([dbc.themes.COSMO]))
dash.register_page(__name__,
                   path='/',
                   title="Homepage",
                   name="Homepage")

layout = dbc.Container(
    html.Div(id = "welcometext",
        children = "Welcome at GLYCOtools.\n"
                   "Navigate on top for the modules you are interested in.\n"
                   "Please note that the dashboard is still under development.",
        style={"whiteSpace" : "pre",
            "font-size" : "26px"}
    )
)
