#!/bin/env python
"""
Custom exceptions.

TODO: How do I turn off printing the last function call in the trace?
Witness this:

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "feature.py", line 186, in __cmp__
    raise InvalidComparison, 'Cannot compare the location of features on opposite strands; they are in different coordinate systems.'
feature.InvalidComparison: Cannot compare the location of features on opposite strands; they are in different coordinate systems.

This is redundant info.  Normal Python exceptions don't print the 'raise'
call.
"""

###############################################################################
# Passing stuff to functions.
###############################################################################

class InvalidArgumentError (StandardError):
    """Indicates an invalid argument to a function."""


###############################################################################
# Other.
###############################################################################

class FeatureNotImplementedError (StandardError):
    """Indicates a feature that has not been implemented."""

