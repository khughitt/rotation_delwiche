"""
Methods for working with HMMER3 CSV output

Keith Hughitt <khughitt@umd.edu>
2012/09/15
"""
import re
import numpy as np
import matplotlib.pyplot as plt

def parse_csv(filepath):
    """Parses a standard hmmersearch CSV table and returns a NumPy recarray"""
    # read first couple lines of file to determine column names and widths
    fp = open(filepath)
    fp.readline()
    nameline = fp.readline()
    dividers = fp.readline()

    # determine column widths and names
    col_widths = []
    col_names = []
    width = 0

    # get a list of the dividers underneath the column names --
    # we will use these to determine the width of each column
    divider_parts = filter(None, re.split("(#?\s*-+\s+)", dividers))

    for col in divider_parts:
        # get column width
        col_size = len(col)

        # add column name and width
        col_names.append(nameline[width: width + col_size].lstrip("#").strip())
        width += col_size
        col_widths.append(col_size)

    # allow for longer descriptions
    col_widths[-1] = 100

    return np.genfromtxt(fp, dtype=None, names=col_names, delimiter=col_widths,
                         autostrip=True)
                         
def plot_evalue_histogram(recarray, filename, title=""):
    """Creates a histogram of E-values for each species"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    n, bins, patches = ax.hist(np.log(recarray['Evalue']), bins=range(-150, 0), 
                               facecolor='limegreen', alpha=0.75)
    plt.xlabel('Log(E-value)')
    plt.ylabel('Frequency')
    title_text = r'$\mathrm{%s\ E-values:}\ n=%d,\ median=%.4g$'
    plt.title(title_text % (title.replace(" ", "\ "), recarray.shape[0], 
                            np.median(recarray['Evalue'])))
    
    plt.savefig(filename + ".png")
