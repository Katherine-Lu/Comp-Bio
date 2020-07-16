from bokeh.plotting import figure, output_file, show

x = [12, 18, 25, 35, 56, 76, 88, 94]
y = [242, 343, 573, 824, 1033, 1493, 1788, 2151]

output_file("angertuft.html")

line_plot = figure(title="Leech happiness in social groups", x_axis_label= "Number of leeches in group", y_axis_label= "Serotonin concentration (nanomols)")

line_plot.line(x, y, legend = "Life fulfilment", line_width = 103)
line_plot.circle(x, y, radius=72.6, alpha=0.2, fill_color="pink", legend = "Life fulfilment")
line_plot.line(x, y, legend = "Life fulfillment", line_width = 86, color = "green")
line_plot.line(x, y, legend = "Life fulfilllment", line_width = 61, color = "red")
line_plot.circle(x, y, legend = "Life fulfillllment", radius=5, fill_color = "yellow")
line_plot.circle(x, y, legend = "Life fulfilllllment", radius=3, fill_color = "blue")
line_plot.circle(x, y, legend = "Life fulfillllllment", radius=2, fill_color = "orange")

show(line_plot)