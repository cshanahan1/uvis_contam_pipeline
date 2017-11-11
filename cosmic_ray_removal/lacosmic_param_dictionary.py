#! /usr/bin/env python

def lacosmic_param_dictionary():
    """Holds the best LACosmic parameters for each filter, primarily
    for white dwarf standard GRW+70. Also works well on other white
    dwarfs and red star P330E.
    If running over a field of stars instead of a single standard, will
    likely need increase 'objlim' param.
    Add filters as needed!

    These are determined using :mod:`run_lacosmic_tester`.

    Parameters:
        nothing

    Returns:
        param_dict : dictionary
            Parameters to be used in LACosmic.
            {'Filter':[sigclip, sigfrac, objlim, niter, sigclip_pf]}

    Outputs:
        nothing
    """
    param_dict = {'F200LP':[10, 0.4, 2, 5, 9.5],
                  'F218W':[7.5, 0.3, 4, 6, 9.5],
                  'F225W':[10, 0.2, 2, 6, 9.5],
                  'F275W':[10, 0.2, 2, 6, 9.5],
                  'F280N':[8, 0.3, 2, 6, 9.5],
                  'F300X':[8, 0.3, 2, 6, 9.5],
                  'F336W':[10, 0.2, 2, 6, 9.5],
                  'F365N':[10, 0.3, 2, 6, 9.5],
                  'F390M':[7, 0.3, 2, 6, 9.5],
                  'F390W':[10, 0.3, 2, 6, 9.5],
                  'F395N':[10, 0.3, 2, 6, 9.5],
                  'F410M':[10, 0.3, 2, 6, 9.5],
                  'F438W':[7.5, 0.2, 2, 6, 9.5],
                  'F343N':[10, 0.3, 2, 6, 9.5],
                  'F373N':[10, 0.3, 2, 6, 9.5],
                  'F467M':[10, 0.3, 2, 6, 9.5],
                  'F469N':[10, 0.3, 2, 6, 9.5],
                  'F475W':[10, 0.3, 2, 6, 9.5],
                  'F502N':[10, 0.3, 2, 6, 9.5],
                  'F547M':[10, 0.3, 2, 6, 9.5],
                  'F555W':[10, 0.3, 2, 6, 9.5],
                  'F606W':[10, 0.2, 2, 6, 9.5],
                  'F631N':[10, 0.3, 2, 6, 9.5],
                  'F645N':[10, 0.3, 2, 6, 9.5],
                  'F656N':[10, 0.3, 2, 6, 9.5],
                  'F657N':[10, 0.3, 2, 6, 9.5],
                  'F658N':[10, 0.3, 2, 6, 9.5],
                  'F665N':[10, 0.3, 2, 6, 9.5],
                  'F673N':[10, 0.3, 2, 6, 9.5],
                  'F680N':[10, 0.3, 2, 6, 9.5],
                  'F689M':[10, 0.3, 2, 6, 9.5],
                  'F763M':[10, 0.3, 2, 6, 9.5],
                  'F775W':[8.5, 0.2, 3, 6, 9.5],
                  'F814W':[8.5, 0.2, 3, 6, 9.5],
                  'F845M':[10, 0.3, 2, 6, 9.5],
                  'F850LP':[10, 0.3, 2, 6, 9.5]}

    return param_dict