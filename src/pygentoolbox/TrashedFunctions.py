def max_of_maximums(dTEmax, dsRNAmax):
    print('Finding max of maximums')
    dmax = {}
    for i in list(dTEmax.keys()):
        for j in list(dsRNAmax.keys()):
            try:
                dTEmax[j]
            except:
                dmax[j] = dsRNAmax[j]
            else:
                dmax[j] = max(dTEmax[j], dsRNAmax[j])

            try:
                dsRNAmax[i]
            except:
                dmax[i] = dTEmax[i]  # if no i in dsRNAmax
            else:  # if scaffold is present in both dictionaries, find the max between them, i and j should be same
                dmax[i] = max(dTEmax[i], dsRNAmax[i])
    for k, v in list(dmax.items()):
        if v == 0:
            print(k)
    return dmax


def max_coords_per_scaffold(d):
    print('Calcualting maximum coordinates per scaffold')
    dmax_coords = {}
    dscaff_max = {}
    for scaff in list(d.keys()):
        # collect all coordinates in list, for each scaffold
        for line in d[scaff]:
            dmax_coords.setdefault(scaff, []).append(int(line.split('\t')[3]))
    for scaff in list(dmax_coords.keys()):
        # get maximum coordinate per scaffold
        dscaff_max[scaff] = max(dmax_coords[scaff])
    return dscaff_max