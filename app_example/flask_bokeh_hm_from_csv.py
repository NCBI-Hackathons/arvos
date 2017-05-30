'''This example demonstrates embedding a standalone Bokeh figure
into a simple Flask application, with a basic HTML web form using Flask.
Presently, this code is a hybrid of examples found at these locations:
    bokeh/examples/embed/simple
    http://bokeh.pydata.org/en/latest/docs/gallery/unemployment.html

For flask to render it, make sure the embed_hm.html file is in the templates folder

Folder structure:
/app_example
    flask_bokeh_heatmap.py
    /templates
        embed_hm.html

To run it type:

    python flask_bokeh_heatmap.py

in the app_example directory, and navigate to:

    http://localhost:5000

to see the heat map
'''
from __future__ import print_function
from math import pi
import pandas as pd
import flask

from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8
from bokeh.models import (
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    BasicTicker,
    PrintfTickFormatter,
    ColorBar,
)
from bokeh.plotting import figure

def get_dataframe_and_axes(fname=None, gene_col_name=None):
    ''' arbitrary data for now '''
    df = pd.read_csv(fname)
    #df.rename(columns=lambda x: x.strip().replace("'",""))
    df['Gene_ID'] = df[gene_col_name].astype(str) 
    df.drop([gene_col_name], axis=1, inplace=True)
    df = df.set_index('Gene_ID')
    df.columns.name = 'Samples' 
    gene_lst = list(df.index)
    sample_lst = list(df.columns)
    df2 = pd.DataFrame(df.stack(), columns=['counts']).reset_index()
    return df2, sample_lst, gene_lst 

def make_heatmap_object(df, sample_lst, gene_lst):
    ''' makes a bokeh figure '''
    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    mapper = LinearColorMapper(palette=colors, low=df['counts'].min(), high=df['counts'].max())

    source = ColumnDataSource(df)

    TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

    # p is a bokeh figure object 
    p = figure(title="Gene Expression",
               x_range=sample_lst, y_range=gene_lst,
               x_axis_location="above", plot_width=600, plot_height=900,
               tools=TOOLS, toolbar_location='below')

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi / 3

    p.rect(x="Samples", y="Gene_ID", width=1, height=1,
           source=source,
           fill_color={'field': 'counts', 'transform': mapper},
           line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                         formatter=PrintfTickFormatter(format="%d"),
                         label_standoff=6, border_line_color=None, location=(0, 0))
    p.add_layout(color_bar, 'right')

    p.select_one(HoverTool).tooltips = [
         ('coord', '@Samples @Gene_ID'),
         ('count', '@counts'),
    ]
    return p


app = flask.Flask(__name__)

@app.route("/")
def home_page():
    """ Simple embedding of a bokeh figure in Flask

    """
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()
    # note that heapmap below is defined under if-name-main block
    script, div = components(heatmap)
    html = flask.render_template(
        'embed_hm2.html',
        plot_script=script,
        plot_div=div,
        js_resources=js_resources,
        css_resources=css_resources
    )
    return encode_utf8(html)

if __name__ == "__main__":
    #print(__doc__)
    filename = 'data/Example_file.csv'
    df, sample_lst, gene_lst = get_dataframe_and_axes(filename, 'Gene ID')
    heatmap = make_heatmap_object(df, sample_lst, gene_lst)
    app.run()

