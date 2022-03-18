def rescale(valuedict,scalemax=256,mincutoff=0):
    """
    Method to rescale a dictionary with numeric values

    :param valuedict: the dictionary with keys and numeric values
    :param scalemax: the max value to scale to
    :param mincutoff: the value below which inputs will be rounded to 0
    :return:
    """

    scalemax = scalemax -1
    smol = min(valuedict.values())
    big = max(valuedict.values())
    for key in valuedict.keys():
        newval = int(valuedict[key]/(big-smol) * scalemax)
        valuedict[key] = newval if newval >= mincutoff else 0
