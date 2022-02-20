##################################################################
############################## 002 ###############################
##################################################################
# H0 fixed
# Data in pH range 0-7.5 to fix X1, KX1
##################################################################

def minner_002(dataframe):
    global H0, X1
    global C_NaOH, K_X1, X2, K_X2

    x_NaOH = dataframe['m']
    data_NaOH = dataframe['H']
    # DEFINE FUNCTION WHICH RETURNS VALUE(S) TO BE MINIMISED, IN THIS CASE SSR 
    def fcn2min(params, x_NaOH, data_NaOH):
        H0 = params['H0']
        C_NaOH = params['C_NaOH']
        X1 = params['X1']
        K_X1 = params['K_X1']
        X2 = params['X2']
        K_X2 = params['K_X2']

        model = ((V0 + Va+ dataframe["m"])*(dataframe["H"]-dataframe["OH"]) 
                 -((V0+Va)*H0)
                 +(dataframe["m"]*C_NaOH) 
                 #- V0*(BT/(1+dataframe["H"]/dataframe["KB"]))
                 - (V0)*(CO2/(1+(dataframe["H"]/(dataframe["K1_LK"]))+dataframe["K2_LK"]/dataframe["H"]))
                 - (V0)*(X1/(1+dataframe["H"]/K_X1))
                 - (V0)*(X2/(1+dataframe["H"]/K_X2))
                 #- (V0)*(X2/(1+dataframe["H"]/K_X3))
                )


        return model - data_NaOH

    # create a set of Parameters
    params = Parameters()
    params.add('H0',   value= H0, vary=False ) #highest [H+] value used as initial estimate
    params.add('C_NaOH',   value= C_NaOH, vary=False ) 
    params.add('X1',   value= X1 )
    params.add('K_X1',   value= K_X1)
    params.add('X2',   value= X2_initial, min = 0)
    params.add('K_X2',   value= K_X2_initial, min = 0)

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x_NaOH, data_NaOH))
    kws  = {'options': {'maxiter':100}}
    result = minner.minimize()
    for _ in range(10):
        minner.minimize()
    # calculate final result
    final = data_NaOH + result.residual

    X1 = result.params.get('X1').value
    X1_std = result.params.get('X1').stderr

    K_X1 = result.params.get('K_X1').value
    K_X1_std = result.params.get('K_X1').stderr

    X2 = result.params.get('X2').value
    X2_std = result.params.get('X2').stderr

    K_X2 = result.params.get('K_X2').value
    K_X2_std = result.params.get('K_X2').stderr
