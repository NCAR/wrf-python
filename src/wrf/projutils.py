from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

from .py3compat import viewitems

def dict_keys_to_upper(d):
    """Return a dictionary with the keys changed to uppercase.
    
    Args:
        
        d (:obj:`dict`): A dictionary.
        
    Returns:
    
        :obj:`dict`: A dictionary with uppercase keys.
    
    """
    return {key.upper() : val for key, val in viewitems(d)}
