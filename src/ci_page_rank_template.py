"""
MIT License

Copyright (c) 2016 Michael Schneider

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import datetime
import os
import sys
import networkx as nx
import InputOutput
import numpy
import argparse
import cPickle

options = {}


def parse_arguments():
    """Specify and parse command line inputs
    """
    global options

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", dest="example", help="An example", required=True)
    options = parser.parse_args()


def build_corroborating_information_graph(data):
    """
    Implement the building of the CI Graph here!

    Parameters
    ----------
    data : Your input data

    Returns
    -------
    """
    pass


def do_page_rank(ci_graph):
    """
    Implement the ci_graph here

    Parameters
    ----------
    ci_graph : CI graph to run pagerank on

    Returns
    -------
    """
    pass


def main():
    """
    Call main program here

    Returns
    -------

    """

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")
