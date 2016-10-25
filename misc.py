#!/bin/env python
"""Dumping ground for stuff that doesn't fit well in other modules.

Some of these should probably be moved to C++ at some point.
"""

def printe (*fields):
    """Shorthand for printing to C{STDERR}; takes same args that you would give to C{print} (comma-delimited objects
    that can be converted to strings).
    """
    print >> sys.stderr, ' '.join (map (str, fields))


