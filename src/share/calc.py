import math
import functions as f
import numpy as np

def number_density(mask,traj,wd):
    print("calculating number densities")
    parameters = {
                  'bins'    : None,
                  'residue' : []
                 }

    center = [float(mask['center'][0]),
              float(mask['center'][1]),
              float(mask['center'][2])
             ]

    binlist = []
    natoms = len(mask['mask'])
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
    steps = mask['volume']/float(parameters['bins'])
    maxbin = int(round(round(steps, 2)))
    bins = {}
    empty = []
    n_density = []

    for i in range(0,maxbin):
        r = (i + 1) * parameters['bins']
        r0 = i * parameters['bins']

        # the last bin can be smaller
        if i + 1 == maxbin:
            if r > float(mask['volume']):
                r = r - (r - mask['volume'])

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
    for frame in traj['frames']:
        bin_tmp = []
        frame = int(frame[0])
        for ai in mask['mask']:
            if mask['mask'][ai][4] == 'HOH' and mask['mask'][ai][2] == 'O':
                coord_index = (frame) * natoms + ai

                for b in bins:
                    if f.euclidian_overlap(center,traj['coords'][coord_index],float(bins[b][1])) == True:
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
        
def EXPfwd(MA1,l1,l2,kT,skip):
    """
        Zwanzig exponential formula
        Need to feed lambdas
    """
    veff1 = 0.0
    veff2 = 0.0
    total = 0.0
    skip = 0  #TEMP
    kT = float(kT)
    for ipt in range(skip,len(MA1) -1):
        for state in range(0,len(l1)):
            veff1 += l1[state] * MA1[ipt][state]
            veff2 += l2[state] * MA1[ipt + 1][state]
          
        dv=veff2-veff1
        veff1=0.0
        veff2=0.0
        total = total + math.exp(-dv/kT)

    avg = total/(len(MA1)-skip)
    dGf = -kT*math.log(avg)
    
    return dGf

#       ! boltzmann's constant beta = 1/rt
#      dgf=0.
#      dgfsum=0.
#      sum=0.
#      veff1=0.
#      veff2=0.
#      dv=0.
#      do ifile=1,nfiles-1
#         do ipt=nskip+1,FEP(ifile)%npts
#            do istate=1,nstates
#               veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#               veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#            end do
#            dv=veff2-veff1
#            veff1=0.
#            veff2=0.
#            sum=sum+exp(-dv/rt)
#         end do
#         sum=sum/real(FEP(ifile)%npts-nskip)
#         dgf()=-rt*dlog(sum)

def EXPbwd(MA1,l1,l2,kT,skip):
    """
        Zwanzig exponential formula
        Need to feed lambdas
    """
    veff1 = 0.0
    veff2 = 0.0
    total = 0.0
    skip = 0  #TEMP
    kT = float(kT)
    MA1 = MA1[::-1]
    for ipt in range(skip,len(MA1) - 1):
        for state in range(0,len(l1)):
            veff1 += l1[state] * MA1[ipt][state]
            veff2 += l2[state] * MA1[ipt + 1][state]
          
        dv=veff1-veff2
        veff1=0.0
        veff2=0.0
        math.exp(-dv/kT)
        total = total + math.exp(-dv/kT)

    avg = total/(len(MA1)-skip)
    dGr = -kT*math.log(avg)
    
    return dGr
#   ! Do the reverse calculation of delta G's according to Zwanzig.
#      sum=0.
#      veff1=0.
#      veff2=0.
#      dv=0.
#      do ifile=nfiles,2,-1
#         do ipt=nskip+1,FEP(ifile)%npts
#            do istate=1,nstates
#               veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#               veff2=veff2+FEP(ifile-1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#            end do
#            dv=veff2-veff1
#            veff1=0.
#            veff2=0.
#            sum=sum+exp(-dv/rt)
#         end do
#         sum=sum/real(FEP(ifile)%npts-nskip)
#         dgr=-rt*dlog(sum)
#      end do

""" ! Overlap Sampling. Lu
     dglu=0.
     dglusum=0.
     sum=0.
     veff1=0.
     veff2=0.
     dv=0.
     do ifile=1,nfiles-1
       ! Note, similar to Zwanzig for the forward calculation
        do ipt=nskip+1,FEP(ifile)%npts
           do istate=1,nstates
              !print *, 'veff1 = FEP(ifile)%lambda(istate) * FEP(ifile)%v(istate,ipt)', FEP(ifile)%lambda(istate), FEP(ifile)%v(istate,ipt), ifile, istate,ipt
              veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
              veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
           end do
           dv=(veff2-veff1)/2
           veff1=0.
           veff2=0.
           sum=sum+exp(-dv/rt)
        end do
        sumf=sum/real(FEP(ifile)%npts-nskip)
        sum=0.

#       ! boltzmann's constant beta = 1/rt
#      dgf=0.
#      dgfsum=0.
#      sum=0.
#      veff1=0.
#      veff2=0.
#      dv=0.
#      do ifile=1,nfiles-1
#         do ipt=nskip+1,FEP(ifile)%npts
#            do istate=1,nstates
#               veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#               veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#            end do
#            dv=veff2-veff1
#            veff1=0.
#            veff2=0.
#            sum=sum+exp(-dv/rt)
#         end do
#         sum=sum/real(FEP(ifile)%npts-nskip)
#         dgf()=-rt*dlog(sum)



        ! Note, similar to Zwanzig for the backward calculation
        do ipt=nskip+1,FEP(ifile+1)%npts
           do istate=1,nstates
              veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
              veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
           end do
           dv=(veff2-veff1)/2
           veff1=0.
           veff2=0.
           sum=sum+exp(dv/rt)
        end do
        sumb=sum/real(FEP(ifile+1)%npts-nskip)

        dglu(ifile)=-rt*dlog(sumf/sumb) ! used later as initial guess for bennetts's constant
        dglusum(ifile+1)=dglusum(ifile)+dglu(ifile)
        sum=0.
     end do

     !print *, 'constant for the bennett formula from overlap sampling', dglu(1) """

def overlap_sampling(MA1,l1,l2,kT,skip):
    """
        Overlap sampling
    """
    veff1 = 0.0
    veff2 = 0.0
    total = 0.0
    skip = 0  #TEMP
    kT = float(kT)
    #MA1 = MA1[::-1]
    for ipt in range(skip,len(MA1) - 1):
        for state in range(0,len(l1)):
            veff1 += l1[state] * MA1[ipt][state]
            veff2 += l2[state] * MA1[ipt + 1][state]
          
        dv=(veff1-veff2)/2
        veff1=0.0
        veff2=0.0
        math.exp(-dv/kT)
        total = total + math.exp(-dv/kT)

    totalf = total / len(MA1)
    veff1 = 0.0
    veff2 = 0.0
    total = 0.0
    skip = 0  #TEMP
    MA1 = MA1[::-1]
    for ipt in range(skip,len(MA1) - 1):
        for state in range(0,len(l1)):
            veff1 += l1[state] * MA1[ipt][state]
            veff2 += l2[state] * MA1[ipt + 1][state]
          
        dv=(veff1-veff2)/2
        veff1=0.0
        veff2=0.0
        math.exp(-dv/kT)
        total = total + math.exp(-dv/kT)

    totalb = total / len(MA1)
    dglu = -kT*math.log(totalf/totalb)

    return dglu

def BAR(MA1,l1,l2,kT,skip,dglu):
    """
        Zwanzig exponential formula
        Need to feed lambdas
    """
#      dgbar=0.
#      dgbarsum=0.
#      sum=0.
#      veff1=0.
#      veff2=0.
#      dv=0.
#      konst=0
    fel=1.0
    felcutoff=0.001
    konst = dglu

    while (fel > felcutoff):
        veff1 = 0.0
        veff2 = 0.0
        total = 0.0
        skip = 0  #TEMP
        kT = float(kT)
        #MA1 = MA1[::-1]
        for ipt in range(skip,len(MA1) - 1):
            for state in range(0,len(l1)):
                veff1 += l1[state] * MA1[ipt][state]
                veff2 += l2[state] * MA1[ipt + 1][state]
            
            dv=veff1-veff2
            veff1=0.0
            veff2=0.0
            math.exp(-dv/kT)
            total = total + 1/(1 + math.exp((dv-konst)/kT))

        totalf = total / len(MA1)

        veff1 = 0.0
        veff2 = 0.0
        total = 0.0
        skip = 0  #TEMP
        MA1 = MA1[::-1]
        for ipt in range(skip,len(MA1) - 1):
            for state in range(0,len(l1)):
                veff1 += l1[state] * MA1[ipt][state]
                veff2 += l2[state] * MA1[ipt + 1][state]
            
            dv=veff1-veff2
            veff1=0.0
            veff2=0.0
            math.exp(-dv/kT)
            total = total + 1/(1 + math.exp((-dv+konst)/kT))

        totalb = total / len(MA1)

#            dgbar(ifile)=-rt*dlog((sumf/sumb)*exp(-konst/rt)*nfnr)
        dgbar = -kT*math.log((totalf/totalb)*math.exp(-konst/kT))
        fel = abs(konst - dgbar)
        konst = dgbar
    
    return dgbar

# !  Bennet's Acceptance Ratio (BAR)
#      dgbar=0.
#      dgbarsum=0.
#      sum=0.
#      veff1=0.
#      veff2=0.
#      dv=0.
#      konst=0
#      fel=1
#      do ifile=1,nfiles-1
#         konst=dglu(ifile)
#         !print *, 'constant for the bennett formula initial', dglu(ifile)
#         felcutoff = 0.001
#         do while (fel > felcutoff)
#            nfnr=real(FEP(ifile)%npts-nskip)/real(FEP(ifile+1)%npts-nskip)
#            nrnf=real(FEP(ifile+1)%npts-nskip)/real(FEP(ifile)%npts-nskip)
#            !print *, 'optimization cutoff begin', 0.001
#            !print *, 'nfnr, nrnf, ifile = ', nfnr, nrnf, ifile

#            do ipt=nskip+1,FEP(ifile)%npts
#               do istate=1,nstates
#                  veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#                  veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile)%v(istate,ipt)
#               end do
#               dv=(veff2-veff1)
#               veff1=0.
#               veff2=0.
#               sum=sum+1/((1+(nfnr*exp((dv-konst)/rt))))
#            end do
#            sumf=sum/real(FEP(ifile)%npts-nskip)

#            sum=0.
#            do ipt=nskip+1,FEP(ifile+1)%npts
#               do istate=1,nstates
#                  veff1=veff1+FEP(ifile)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
#                  veff2=veff2+FEP(ifile+1)%lambda(istate)*FEP(ifile+1)%v(istate,ipt)
#               end do
#               dv=(veff2-veff1)
#               veff1=0.
#               veff2=0.
#               sum=sum+1/((1+(nrnf*exp((-dv+konst)/rt))))
#            end do
#            sumb=sum/real(FEP(ifile+1)%npts-nskip)


#            dgbar(ifile)=-rt*dlog((sumf/sumb)*exp(-konst/rt)*nfnr)
#            sum=0.
#            fel=abs(konst-dgbar(ifile))
#            konst=dgbar(ifile)
#            !print *, 'optimization cutoff end', ifile, fel
#         end do

#         !print *, 'constant for the bennett formula end', dglu(ifile)
#         dgbarsum(ifile+1)=dgbarsum(ifile)+dgbar(ifile)
#         !print *, 'The optimization constant is', fel
#         fel=1

#      end do