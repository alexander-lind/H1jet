#!/usr/bin/env python3 

# Usage: 
#    ./h1jet | python PlotH1jet.py 
# or: 
#    ./h1jet -o result.out 
#    python PlotH1jet.py result.out 

import sys 
import select
import matplotlib
from matplotlib import pyplot as plt

# Containers 
bin_min    = []
bin_med    = []
bin_max    = []
dsigma_dpt = []
sigmapt    = []
collider   = "pp"
roots      = 0
process    = r"$pp \rightarrow H + $jet"
model      = "SM"
pdf_name   = ""
log        = False 
legend_x0  = 0.55 
pdf_line   = False
hist_line  = False 

def parseLine(line):
    global pdf_line
    global hist_line
    rline = line.rstrip() 
    if (pdf_line == True):
        pdf_name = rline.split()[0] 
        pdf_line = False
    elif rline.startswith(' collider       = ppbar'):
        collider = r"$p\bar{p}$"
    elif rline.startswith(' process        = Z'):
        process = r"$pp \rightarrow Z + $jet"
    elif rline.startswith(' process        = bbH'):
        process = r"$b\bar{b} \rightarrow H + $jet"
    elif rline.startswith(' process        = user'):
        process = "User process" # Change this to your process  
        model = "User model" # Change this to your model 
        legend_x0 = 0.40 
    elif rline.startswith(' model          = SUSY'):
        model = "SUSY"
        legend_x0 = 0.50 
    elif rline.startswith(' model          = SM with top partners'):
        model = "SM + top partner"
        legend_x0 = 0.30 
    elif rline.startswith('LHAPDF'):
        pdf_line = True
    elif rline.startswith(' log            = T'):
        log = True
    elif rline.startswith(' roots(GeV)'):
        roots = int(float(rline.split()[2])/1000.0)
    elif rline.startswith(' # cols'):
        hist_line = True
    elif rline.startswith('   ') and (hist_line == True): 
        data = [float(s) for s in rline.split()]
        bin_min.append(data[0])
        bin_med.append(data[1])
        bin_max.append(data[2])
        dsigma_dpt.append(data[3])
        sigmapt.append(data[4])

def main(): 

    if (len(sys.argv) > 1):
        # Read from input file
        with open(sys.argv[1]) as fp: 
            line = fp.readline() 
            while line: 
                parseLine(line) 
                line = fp.readline()
    elif (select.select([sys.stdin,],[],[],0.0)[0]):
        # Read from stdin 
        for line in sys.stdin: 
            parseLine(line)
    else: 
        print("No input provided")
        sys.exit()

    # Handle Plotting 
    font_size = 11

    bins = len(bin_min)

    x_min = bin_min[0]
    x_max = bin_max[bins - 1]

    # Set Up Figure 
    fig, ax = plt.subplots(num = 1, figsize = (5, 4), dpi = 100)

    # Use scientific notation in axes ticks 
    ax.ticklabel_format(axis = "y", style = "sci", scilimits = (0, 0)) 

    # Show Minor Ticks 
    plt.minorticks_on()

    # Set ticks on the inside, and on all sides 
    ax.tick_params(axis = "y", direction = "in", which = "both", bottom = True, top = True, left = True, right = True)
    ax.tick_params(axis = "x", direction = "in", which = "both", bottom = True, top = True, left = True, right = True)

    # Plot data 
    plt.plot(bin_med, dsigma_dpt, color = "blue", linewidth = 1, linestyle = 'solid')

    # Logarithmic y scale 
    plt.yscale("log")
    #plt.xscale("log")

    # Show all powers of 10 on the y axis 
    locmaj = matplotlib.ticker.LogLocator(base = 10, numticks = 15) 
    ax.yaxis.set_major_locator(locmaj)

    # Show minor ticklabels 
    locmin = matplotlib.ticker.LogLocator(base = 10.0, subs = (0.2, 0.4, 0.6, 0.8), numticks = 15)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    # Labels 
    ax.set_ylabel(r"$\mathrm{d}\sigma / \mathrm{d}p_{T}$ [fb/GeV]", fontsize = font_size)
    ax.set_xlabel(r"$p_{T}$ [GeV]", fontsize = font_size)
    if (log): 
        ax.set_xlabel(r"$\ln\ p_{T}$ [GeV]", fontsize = font_size)

    # Positioning of Labels 
    ax.yaxis.set_label_coords(-0.110, 0.815)
    ax.xaxis.set_label_coords(0.920, -0.075)

    # Range of x-axis 
    plt.xlim(x_min, x_max)

    # Text on Plot 
    plt.text(legend_x0, 0.90, r"H1jet ", transform = ax.transAxes, fontsize = font_size)
    plt.text(legend_x0, 0.82, process, transform = ax.transAxes, fontsize = font_size)
    plt.text(legend_x0, 0.74, model + ", " + collider + r"$, \sqrt{s} = $" + str(roots) + " TeV", transform = ax.transAxes, fontsize = font_size)
    plt.text(legend_x0, 0.66, pdf_name, transform = ax.transAxes, fontsize = font_size)

    plt.tight_layout()

    # Show Plot 
    plt.show()

    # Save Plot 
    fig.savefig("graph.eps")
    
    print("Saved figure in graph.eps")

if __name__ == "__main__": 
    main() 

