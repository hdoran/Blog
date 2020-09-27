import plotly.express as px
import numpy as np
p = np.arange(0,1,.0001)	 
result = 400*p**4
fig =px.scatter(x=1-p, y=result)
fig.update_layout(title='Brute Force Search Over Various Values of P%',
	xaxis_title='Value for P',
	yaxis_title='Ball Location at L_4')
fig.add_shape(
	# Line Horizontal
		type="line",
		x0=0,
		y0=2.125,
		x1=1,
		y1=2.125,
		line=dict(
			color="Green",
			width=4,
		),
)
fig.update_shapes(dict(xref='x', yref='y'))		
fig.write_html("C:/GitStuff/golf_plot.html")	
fig.show() 