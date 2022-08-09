import pandas as pd
import plotly.express as px

root = "./endpoints/"

final_df = pd.read_csv(root + "ScoreDF" +".csv")

fig = px.scatter(final_df, "Position", "Scores", color = "Scores", hover_data = ["Residue"])

fig.show()