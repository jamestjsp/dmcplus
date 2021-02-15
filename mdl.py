def get_dmc_model(file):
    '''
    A function to parse Apetch DMCPlus model
    '''
    # mdl file coments
    sComment = next(file).strip()
    #Number of files and the file names
    sBuf = int(next(file).strip())
    for _ in range(sBuf):
        next(file)
    # Number of independents
    NumberOfIndependents = int(next(file).strip())
    # Number of dependents
    NumberOfDependents = int(next(file).strip())
    # Number of coefficients
    NumberOfCoefficients = int(next(file).strip())
    # Time to steady-state
    next(file) #Skip one line
    SteadyStateTime = float(next(file).strip())
    # DMCplus model style flag
    sBuf = int(float(next(file).strip()))
    bNewStyle = True if sBuf == 9896 else False
    # read tagnames
    Independents = list()
    for _ in range(NumberOfIndependents):
        sBuf = next(file).strip().split()[1]
        Independents.append(sBuf)
    Dependents = list()
    for _ in range(NumberOfDependents):
        sBuf = next(file).strip().split()[1]
        Dependents.append(sBuf)
    isRamp = dict()
    dGain = dict()
    fir_curves = dict()
    for dep in range(NumberOfDependents):
        dGain[Dependents[dep]] = {}
        fir_curves[Dependents[dep]] = {}

        # Read & thrash the next dependent variable header: save the ramp status
        isRamp[Dependents[dep]] = True if int(next(file).strip().split()[1]) else False
        for _ in range(10):
            next(file)
        for ind in range(NumberOfIndependents):
            Coefficient = list()
            # Read & Store the next independent variable curve
            # Tagname, Eng Units, and double precision SS gain
            dGain[Dependents[dep]][Independents[ind]] = float(next(file).strip().split()[-1])
            # Model coefficients
            while len(Coefficient) < NumberOfCoefficients:
                Coefficient.extend([float(i) for i in next(file).strip().split()])
            fir_curves[Dependents[dep]][Independents[ind]] = Coefficient
            continue
        continue
    DMCPlusModel = {'Coment':sComment, 
                    'NumberOfIndependents':NumberOfIndependents, 
                    'NumberOfDependents':NumberOfDependents,
                    'NumberOfCoefficients':NumberOfCoefficients,
                    'SteadyStateTime':SteadyStateTime,
                    'NewStyle':bNewStyle,
                    'Independents':Independents,
                    'Dependents':Dependents,
                    'isRamp':isRamp,
                    'dGain':dGain,
                    'Coefficients':fir_curves
                    }
    return DMCPlusModel

if __name__ == '__main__':
    import  json
    infile = 'mdl/2x2_model.mdl'
   
    with open(infile, 'r') as f:
        curves = get_dmc_model(f)
    
    outfile = 'json' + infile.rsplit('.',1)[0][3:] + '.json'
    
    with open(outfile, "w") as outfile:
        json.dump(curves, outfile)


