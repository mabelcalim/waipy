from __future__ import absolute_import
import os


def example_data_path(filename=None):
    """Return the full path to the example data directory."""
    if filename is None:
        filename = ''
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        'example_data', filename)


__all__ = ['example_data_path']

