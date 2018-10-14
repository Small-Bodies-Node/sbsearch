# Licensed with the 3-clause BSD license.  See LICENSE for details.
"""utility closet"""
import numpy as np
import astropy.units as u


def assemble_sql(cmd, parameters, constraints):
    """Assemble a SQL statement.

    Parameters
    ----------
    cmd : string
        Left-hand side of the statement.

    parameters : list
        Parameters for substitution (via sqlite3 parameter
        substitution).

    constraints : list of tuple
        Each constraint is a SQL expression and an optional parameter
        for subsitution into the expression.  If no parameter is used,
        set it to ``None``.

    """
    if len(contraints) > 0:
        expr, params = list(zip(constraints))
        cmd = cmd + ' WHERE ' + ' AND '.join(expr)
        parameters.extend(params)
    return cmd, parameters


def date_constraints(jd_start, jd_stop):
    """Add date constraints for assemble_sql()."""
    constraints = []
    if jd_start is not None:
        constraints.append(('jd>=?', jd_start))

    if jd_stop is not None:
        constraints.append(('jd<=?', jd_stop))

    return constraints


def iterate_over(cursor):
    while True:
        rows = cursor.fetchmany()
        if not rows:
            return
        for row in rows:
            yield row
