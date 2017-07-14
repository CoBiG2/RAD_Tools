#!/usr/bin/python3
# Copyright 2017 Diogo N Silva <o.diogosilva@gmail.com>
# plot_structure is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# plot_structure is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with plot_structure. If not, see <http://www.gnu.org/licenses/>.

# plot_structure uses plotly library to produce Structure interactive Plots


# Plotly imports
from plotly.offline import plot
import plotly.graph_objs as go
from plotly import tools

# Other imports
import argparse
import numpy as np
import os
import operator
import colorlover as cl

parser = argparse.ArgumentParser(description="Produces interactive structure"
                                             " plots.")
parser.add_argument("-in", dest="infile", nargs="*", help="meanQ file")
parser.add_argument("-o", dest="output_file", help="Name of plot file")
parser.add_argument("-p", dest="popfile", help="")

arg = parser.parse_args()


def parse_q(qfile):
    """
    Parses a meanQ file from fastStructure and returns a list of lists, with
    the values of each taxon in each element of the root list
    """

    fh = open(qfile)

    return np.array([[float(i) for i in x.strip().split()] for x in
                      fh.readlines()])


def parse_pops(popfile):
    """
    Parses a population file and returns a list in which the index of the taxon
    names correspond to their positions in the qfile
    """

    fh = open(popfile)

    return np.array([[x.strip()] for x in fh.readlines()])


def html_creator(plot_div, output_file):

    fh = open(output_file + ".html", "w")

    template = """
<html><head><meta charset="utf-8" /></head><body>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<div id="htmlwidget_container">
<input id="inputText" value></input>
<button id="buttonSearch">
Search
</button>
<script>document.getElementById("buttonSearch").addEventListener("click", function()
  {
    var i = 0;
    var j = 0;
    var found = [];
    var myDiv = document.getElementsByClassName("xtick")
    for (i = 0; i < myDiv.length; i++) {
    	myDiv[i].style.fontWeight = "normal";
    	myDiv[i].childNodes[0].style.fill = "black";

    	if (document.getElementById("inputText").value !== "" && myDiv[i].textContent.indexOf(document.getElementById("inputText").value) !== -1) {
    		myDiv[i].style.fontWeight = "bold";
    		myDiv[i].childNodes[0].style.fill = "red";
    	}
    }
    Plotly.Fx.hover(myDiv, found);
  }
);</script>
</div>
%s
</body></html>
    """ % plot_div

    fh.write(template)


def order_qvals(qvals, pops):

    ranking = []

    for p, q in enumerate(qvals):
        ranking.append([list(q).index(max(q)), float(max(q)), str(pops[p][0])])

    ranking2 = [x[2] for x in reversed(sorted(ranking))]

    return ([[i for i in pops].index([x]) for x in ranking2],
            ranking2,
            sorted(ranking))

def plot_q(qvals, output_name, pops=None):
    """
    Plots the qvals. If pops are provided, plot them as well
    """

    # Get number of structure Plots
    nplots = len(qvals)

    c = None

    fig = tools.make_subplots(rows=nplots, cols=1,
                              shared_xaxes=True,
                              subplot_titles=[os.path.basename(x) for x in qvals],
                              vertical_spacing = 0.05)


    for j, (name, qval) in enumerate(qvals.items()):

        # Get K number
        k = qval.shape[1]

        if not c:
            c = cl.scales[str(k)]["qual"]["Set2"]

        # If pops have not been provided, produce empty x axis labels
        if not any(pops):
            x_labels = np.array([[str(x)] for x in range(qval.shape[0])])
        else:
            x_labels = pops

        if j == 0:
            rank, pops, assignments = order_qvals(qval, pops)

            with open("group_assignment.csv", "w") as fh:
                for k in assignments:
                    fh.write("{}, Cluster_{}\n".format(k[2], k[0]))

        data = []

        for p, i in enumerate(qval.T):

            current_bar = go.Bar(
                x = pops,
                y = [i[x] for x in rank],
                #y = i,
                legendgroup = "group{}".format(p),
                text = ["Assignment: {}%".format(x * 100) for x in i],
                marker=dict(
                    color=c[p],
                    line=dict(
                        color='rgb(8,48,107)',
                        width=0.9,
                    ),
                ),
                name = "K {}".format(p),
                showlegend = True if j == (0) else False
            )

            fig.append_trace(current_bar, j + 1, 1)

    fig["layout"].update(barmode="stack", bargap=0)
    #plot(fig, filename=output_name)
    p_html = plot(fig, include_plotlyjs=False, output_type='div')

    html_creator(p_html, output_name)


def main():
    qfiles = arg.infile
    outfile = arg.output_file
    popfile = arg.popfile

    qvals = {}
    for qf in reversed(qfiles):
        qvals[os.path.basename(qf)] = parse_q(qf)

    pop_vals = parse_pops(popfile)

    plot_q(qvals, outfile, pop_vals)

main()
