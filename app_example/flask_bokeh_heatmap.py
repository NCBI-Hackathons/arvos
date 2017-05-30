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


def get_dataframe_and_axes():
    ''' arbitrary data for now '''
    from bokeh.sampledata.unemployment1948 import data
    data['Year'] = data['Year'].astype(str)
    data = data.set_index('Year')
    data.drop('Annual', axis=1, inplace=True)
    data.columns.name = 'Month'
    years = list(data.index)
    months = list(data.columns)
    # reshape to 1D array or rates with a month and year for each row.
    df = pd.DataFrame(data.stack(), columns=['rate']).reset_index()
    return df, years, months

def make_heatmap_object(df, years, months):
    ''' makes a bokeh figure '''
    # this is the colormap from the original NYTimes plot
    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    mapper = LinearColorMapper(palette=colors, low=df.rate.min(), high=df.rate.max())

    source = ColumnDataSource(df)

    TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"

    # p is a bokeh figure object 
    p = figure(title="US Unemployment ({0} - {1})".format(years[0], years[-1]),
               x_range=years, y_range=list(reversed(months)),
               x_axis_location="above", plot_width=900, plot_height=400,
               tools=TOOLS, toolbar_location='below')

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = pi / 3

    p.rect(x="Year", y="Month", width=1, height=1,
           source=source,
           fill_color={'field': 'rate', 'transform': mapper},
           line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                         formatter=PrintfTickFormatter(format="%d%%"),
                         label_standoff=6, border_line_color=None, location=(0, 0))
    p.add_layout(color_bar, 'right')

    p.select_one(HoverTool).tooltips = [
         ('date', '@Month @Year'),
         ('rate', '@rate%'),
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
        'embed_hm.html',
        plot_script=script,
        plot_div=div,
        js_resources=js_resources,
        css_resources=css_resources
    )
    return encode_utf8(html)

if __name__ == "__main__":
    print(__doc__)
    df, years, months = get_dataframe_and_axes()
    heatmap = make_heatmap_object(df, years, months)
    app.run()

