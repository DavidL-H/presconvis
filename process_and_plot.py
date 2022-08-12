import pandas as pd
import plotly.express as px

root = "./endpoints/"

final_df = pd.read_csv(root + "UniRef50_Q8A7T2_final_df" +".csv")

fig = px.scatter(final_df, "Position", "Scores", color = "Scores", hover_data = ["Residue"])
fig.show()
fig.write_html("./endpoints/plot.html")