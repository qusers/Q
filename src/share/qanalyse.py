def number_density():
    print("calculating number densities")
    parameters = {
                  'bins'    : None,
                  'residue' : []
                 }

    center = [float(Trajectory.center[0]),
              float(Trajectory.center[1]),
              float(Trajectory.center[2])
             ]

    binlist = []

    # read in parameters from file
    # probably needs some dir logic
    with open('number_density.calc') as infile:
        for line in infile:
            line = line.split()
            if len(line) < 1:
                continue
            if not line[0] in parameters:
                print("FATAL: parameter {} not found in inputfile".format(line[0]))
                sys.exit()

            else:
                if type(parameters[line[0]]) == list:
                    parameters[line[0]].append(line[1])

                else:
                    parameters[line[0]] = int(line[1])

    # Now we need calculate whether the atom falls within a bin.
    # First construct the bins and calculate the volumes
    steps = Trajectory.volume/float(parameters['bins'])
    maxbin = int(round(round(steps, 2)))
    bins = {}
    empty = []
    n_density = []

    for i in range(0,maxbin):
        r = (i + 1) * parameters['bins']
        r0 = i * parameters['bins']

        # the last bin can be smaller
        if i + 1 == maxbin:
            if r > float(Trajectory.volume):
                r = r - (r -Trajectory.volume)

        V = (4 * math.pi * (r ** 3))/3
        bins[i] = [r0,r,V,[]]
        binlist.append(i)

    V_tmp = {}
    for b in bins:
        bins_tmp = bins
        tmp = bins[b][2]
        for i in binlist[0:b]:
            tmp = tmp - bins[i][2]
        V_tmp[b] = tmp

    for tmp in V_tmp:
        bins[tmp][2] = V_tmp[tmp]

    # Now loop over the frames
    for frame in Trajectory.frames:
        bin_tmp = []
        frame = int(frame[0])
        for ai in Trajectory.maskPDB:
            if Trajectory.maskPDB[ai][4] == 'HOH' and Trajectory.maskPDB[ai][2] == 'O':
                coord_index = (frame) * int(Trajectory.natoms) + ai

                for b in bins:
                    if f.euclidian_overlap(center,Trajectory.coords[coord_index],float(bins[b][1])) == True:
                        bins[b][3].append(1)

                    else:
                        bins[b][3].append(0)

        # Calculate number of atoms
        for b in binlist:
            nats = np.sum(bins[b][3])
            for i in range(0,b):
                nats = nats - np.sum(bins[i][3])

            # Calculate the density
            Rho = nats/bins[b][2]
            bin_tmp.append(Rho)


        n_density.append(bin_tmp)
        # Reset the list
        for b in binlist:
            bins[b][3] = []

    data = np.asarray(n_density)
    avg_data = np.mean(data,0)
    sdv_data = np.std(data,0)

    for i in range(0, len(avg_data)):
        print('bin:   {}   avg: {:.4f} +/- {:.4f}'.format(i,
                                                         avg_data[i],
                                                         sdv_data[i]))
