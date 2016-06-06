"""
Some custom exceptions for the NEUS application.
"""

class SupportError(Exception):
    """
    Implements a custom exception that indicates that a call to Partition yielded a coorinate with zero support on that
    """
    def __init__(self, value):
        """
        Construct instance of SupportError.
        """
        self.value = value
    def __str__(self):
        """
        Return string representation of error message. 
        """
        return repr(self.value)
