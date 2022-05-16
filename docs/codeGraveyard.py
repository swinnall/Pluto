# general function for isotherm and surface excess modules
#def getFile(fileDIR,equipParams):

#    with open(fileDIR, newline = '') as f:
#        reader = csv.reader(f, delimiter=equipParams[2])
#        data = list(reader)

    # filter out empty lines (seemingly randomly introduced in Nima files)
#    data = [x for x in data if x != []]

#    return data

# just for getting BeagleHole ellipsometry files, each row needs to be split so don't integrate with getFile
def getEllipsometryFile(fileDIR):

    df = []
    with open(fileDIR, newline = '') as f:
        rdr = csv.DictReader(filter(lambda row: row[0]!='#', f))
        for row in rdr:
            df.append( row[None][0].split())

    return df


def ellipsModel(AOI, par):

    # convert AOI to radians
    theta = []
    for i in range(len(AOI)):
        theta.append(AOI[i] * np.pi / 180)


    # wavelength of laser
    lmbda = 632.8E-9


    # refractive index of air
    N0 = 1


    # film 1
    k1 = 0
    n1 = par[0]
    d1 = par[1]

    Nt1 = n1 - 1j*k1


    # film 2
    k2 = 0
    n2 = par[2]
    d2 = par[3]

    Nt2 = n2 - 1j*k2


    # substrate
    k3 = 0
    n3 = 1.337

    Nt3 = n3 - 1j*k3


    # unknown
    ct2 = -12.5663706143592*1j


    Rho_model = []
    for i in range(len(theta)):

        C0    = np.cos(theta[i])
        S0    = np.sin(theta[i])

        C1    = np.sqrt( 1 - (np.sin(theta[i])*(N0/Nt1))**2  )
        rp10  = (Nt1*C0 - N0*C1) / (Nt1*C0 + N0*C1)
        rn10  = (N0*C0 - Nt1*C1) / (N0*C0 + Nt1*C1)

        C2    = np.sqrt( 1 - (np.sin(theta[i])*(N0/Nt2))**2  )
        rp21  = (Nt2*C1 - Nt1*C2) / (Nt2*C1 + Nt1*C2)
        rn21  = (Nt1*C1 - Nt2*C2) / (Nt1*C1 + Nt2*C2)

        T2    = np.exp(ct2*Nt2*C2*d2/lmbda)

        C3    = np.sqrt( 1 - (np.sin(theta[i])*(N0/Nt3))**2  )
        rp32  = (Nt3*C2 - Nt2*C3) / (Nt3*C2 + Nt2*C3)
        rn32  = (Nt2*C2 - Nt3*C3) / (Nt2*C2 + Nt3*C3)

        crp21 = (rp21 + rp32*T2) / (1 + rp21*rp32*T2)
        crn21 = (rn21 + rn32*T2) / (1 + rn21*rn32*T2)

        T1    = np.exp(ct2*Nt1*C1*d1/lmbda);
        crp10 = (rp10 + crp21*T1) / (1 + rp10*crp21*T1)
        crn10 = (rn10 + crn21*T1) / (1 + rn10*crn21*T1)


        if math.isnan(crp10/crn10) == True:
            pass
        else:
            Rho_model.append((crp10/crn10).real)


    Psi = []
    for i in range(len(Rho_model)):
        Psi.append( np.arctan(abs(Rho_model[i])) )


    Delta = []
    for i in range(len(Rho_model)):
        Delta.append( np.angle(Rho_model[i]) * 180/np.pi )

        if Delta[i] < 0:
            Delta[i] == Delta[i] + 360

    return Rho_model




# least square condition
def leastsquare(Rho_exp, Rho_model):

    diffSq = []
    for i, j in zip(Rho_exp,Rho_model):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)



# residual function for genetic algorithm
def residuals(par, AOI, Rho_exp):

    # where ellipsModel is the ellipsometry theory function
    Rho_model = ellipsModel(AOI, par)

    return leastsquare(Rho_exp, Rho_model)



def parameterOptimisation(AOI, psi_AOI, delta_AOI):

    # update variable names for easier processing
    AOI       = AOI.get(0)
    Psi_exp   = psi_AOI.get(0)
    Delta_exp = delta_AOI.get(0)

    # calculate experimental rho
    Rho_exp = []
    for i in range(len(Psi_exp)):
        Rho_exp.append( np.tan(Psi_exp[i]) * cm.exp( 1j * Delta_exp[i] ).real )

    # define input parameters; [n1, d1, n2, d2] where x1 are film and x2 are water
    par = [1, 1, 1.33, 0]

    # associated parameter bounds; could fix pars by defining in model function
    bounds = [(0,2), (0,2), (1.33,1.33), (0.0,0.0)]

    # genetic algorithm
    geneticOutput = opt.differential_evolution(residuals, bounds, args=(AOI, Rho_exp), maxiter=1000) # might need args=*par or make par global

    # parameter solution
    solution = geneticOutput.x
    print("\nParameter solution (n1, d1, n2, d2):\n %s" %solution)

    # associated cost
    lstsq = residuals(solution, AOI, Rho_exp)
    print("\nCost of chosen solution: %f" %lstsq)

    # number of iterations
    print("\nNumber of iterations: %s" %geneticOutput.nit)

    # bool of success
    print("\nOptimisation status: %s" %geneticOutput.success)

    # termination message
    print("\nTermination message: %s" %geneticOutput.message)

    return



A, lmbda1, phi1, B, lmbda2, phi2, c

## Parameters of petri dish test
# A      = 0.2
# Lmbda1 = 6475
# phi1   = -0.972
# B      = 0.1040
# lmbda2 = 79744
# phi2   = 3.140
# c      = 180.23


" Fit Sine curve to sinusoidal ellipsometry data "

#if analysisType == 'parOpt':
#    parameterOptimisation(AOI, psi_AOI, delta_AOI)

if analysisType == 'fitSin':
    fitSine(t, delta_t)

def sinModel(x,par):

    # initial guess parameters; amplitude, wavenumber, phase
    A      = par[0]
    lmbda1 = par[1]
    phi1   = par[2]
    B      = par[3]
    lmbda2 = par[4]
    phi2   = par[5]
    c      = par[6]

    # basic model, might add second sine function and/or time dependence
    sinY_Model = A*np.sin(2*np.pi*np.array(x)/lmbda1 + phi1) + B*np.sin(2*np.pi*np.array(x)/lmbda2 + phi2) + c

    return sinY_Model


# least square condition
def leastsquare(delta_exp, sinY_Model):

    diffSq = []
    for i, j in zip(delta_exp,sinY_Model):
        diffSq.append( np.power(i-j,2) )

    return np.sum(diffSq)


# residual function for genetic algorithm
def residuals(par, t, delta_exp):

    # where sinModel is the sinusoidal function
    sinY_Model = sinModel(t, par)

    return leastsquare(delta_exp, sinY_Model)


def fitSine(t, delta_t):

    # update variable names for easier processing
    t         = t.get(0)
    Delta_exp = delta_t.get(0)


    # high pass filter: remove values less than 0.01 (final smoothing of function)
    tolerance  = 0.01
    smoothT = []; smoothList = []
    for i in range(0,len(Delta_exp)-1):

        if abs(Delta_exp[i] - Delta_exp[i+1]) < tolerance and Delta_exp[i] < 180.60 and Delta_exp[i] > 179.81:
            smoothList.append(Delta_exp[i])
            smoothT.append(t[i])

    # define boundaries of k; 2pi/lambda
    lmbda1_lower = 5000
    lmbda1_upper = 15000
    lmbda1_guess = 6.429e+03
    lmbda2_lower = 60000
    lmbda2_upper = 80000
    lmbda2_guess = 7.98488484e+04

    # define input parameters; [A, lmbda1, phi1, B, lmbda2, phi2, c]
    #par = [0.4, lmbda1_guess, 0, 180]
    par = [0.4, lmbda1_guess, 0, 0.1, lmbda2_guess, 0, 180]

    # associated parameter bounds; could fix pars by defining in model function
    #bounds = [(0.2,0.5), (lmbda1_lower,lmbda1_upper), (-np.pi,np.pi), (179.5,180.5)]
    bounds = [(0.2,0.5), (lmbda1_lower,lmbda1_upper), (-np.pi,np.pi), (0.0,0.5), (lmbda2_lower,lmbda2_upper), (-np.pi,np.pi), (178,182)]

    # genetic algorithm
    geneticOutput = opt.differential_evolution(residuals, bounds, args=(smoothT, smoothList), maxiter=1000) # might need args=*par or make par global

    # parameter solution
    solution = geneticOutput.x
    print("\nParameter solution (A, lmbda1, phi1, B, lmbda2, phi2, c):\n %s" %solution)

    # associated cost
    lstsq = residuals(solution, t, Delta_exp)
    print("\nCost of chosen solution: %f" %lstsq)

    # number of iterations
    print("\nNumber of iterations: %s" %geneticOutput.nit)

    # bool of success
    print("\nOptimisation status: %s" %geneticOutput.success)

    # termination message
    print("\nTermination message: %s" %geneticOutput.message)


    # generate second sine set
    sinY_Model = sinModel(smoothT,solution)

    # generate subtracted data
    DeltaSmoothSubtracted = []
    for i in range(len(smoothT)):
        DeltaSmoothSubtracted.append( smoothList[i] - sinY_Model[i] + solution[6])

    plt.plot(t, Delta_exp, '.', label='Raw Data', c='#1643A2')
    plt.plot(smoothT, smoothList, '-', label='Smoothed Data', c='#3CA8AB')
    plt.plot(smoothT, sinY_Model, label='Fit', linewidth=3, c='#E9BC5D')
    #plt.plot(smoothT, DeltaSmoothSubtracted, label='Subtracted', linewidth=1.5, c='red') # , c='#E9BC5D'
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Delta (deg)')
    plt.show()


    return


## this is to calibrate the periodicity subtraction of the data
## Create simpleFit option that analyses the periodicity of sample data - could extend to comparison between multiple fits
## this should then generate parameters
## generate sinusoidal data (general function) with these parameters
## plot original data, fit, and subtracted data
## ask user if they would like to update default parameters with the new ones (in practice it's most like the program's default will be set and rarely updated as the same should be used continuously)
## ask user if they want to subtract the sin curve from this data

### Points to think about
## maybe each time the experimental data is read it should be analysed for periodicity [have to select region]
## would have to have a noise filter
## could callibrate against Andreas' data if he still has it
