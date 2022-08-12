'''
Shiny python app for the PResConVis tool
David Lennox-Hvenekilde
220812
'''
from shiny import App, render, ui
from uniprotAPI import all_fasta_query, uniprot_query, all_df_query
from clustalomegaAPI import post_alignment_request, clustalo_alignment, root
from aa_variability_msa import clustalo_to_matrix, clustalo_df_formating, final_MSA_df
import pandas as pd
import time
import os
import pandas as pd
import plotly.express as px


app_ui = ui.page_fluid(
    ui.panel_title("PResConVis"),
    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.input_slider("n", "N", 0, 100, 20),
            ui.output_text_verbatim("txt"),
            ui.markdown("""
        **Work in progress**  
        Quick and simple PRotein RESidue CONservation and VISualization tool.

        <a href="https://github.com/DavidL-H/presconvis/pulse" alt="Activity"><img src="https://img.shields.io/github/commit-activity/m/DavidL-H/presconvis" /></a>
        ![Python 3](https://img.shields.io/badge/Language-Python_3-red.svg)
        ![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)

        This tool takes a single protein sequence and user input, grabbing n proteins from uniprot (UnirefXX), performing MSA, and visualizing both a tree for total conservation, but also a 2d sequence plot, showing indivudual residue conservation.

        The front end will be build in with the newly launched Shiny for Python.
        Persistent data storage wil MySQL

        This tools accesses protein sequences and their meta data from the **Uniprot REST API**:  
        https://rest.uniprot.org/  
        https://www.uniprot.org/help/api  

        Multiple sequence alignments are run remotely with Clustal Omega via the **EMBL-EBI REST API**:  
        https://www.ebi.ac.uk/Tools/services/rest/clustalo  
        https://www.ebi.ac.uk/Tools/common/tools/help/  

            """),
        ),
        ui.panel_main(
            #ui.output_plot("viz"),
            ui.output_table("table_data"),
            #ui.include_HTML(root + "plot.html")
      ),
    )
    
)


def server(input, output, session):
    # Renders amino acid overview table
    @output
    @render.table
    # Replace with reactive data frame based on user input
    def table_data():
        df = pd.read_csv(root + "UniRef50_Q8A7T2_final_df" +".csv", index_col=0)
        return df
    
    @output
    @render.plot
    def viz():
        df = pd.read_csv(root + "UniRef50_Q8A7T2_final_df" +".csv", index_col=0)
        fig = px.scatter(df, "Position", "Scores", color = "Scores", hover_data = ["Residue"])
        return fig

    @output
    @render.text
    def txt():
        return f"n*2 is {input.n() * 2}"
    
app = App(app_ui, server)
