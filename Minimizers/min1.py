
##################################################################
############################## 001 ###############################
##################################################################
# Data in pH range 0-5 used to constrain H0, X1, K_X1
##################################################################

def minner_001(dataframe):
    global H0, X1, K_X1

    x_NaOH = dataframe['m']
    data_NaOH = dataframe['H']
    # DEFINE FUNCTION WHICH RETURNS VALUE(S) TO BE MINIMISED, IN THIS CASE SSR 

    def fcn2min(params, x_NaOH, data_NaOH):
        H0 = params['H0']
        C_NaOH = params['C_NaOH']
        X1 = params['X1']
        K_X1 = params['K_X1']

        model = ((V0 + Va+ dataframe["m"])*(dataframe["H"]-dataframe["OH"]) 
                 -((V0+Va)*H0)
                 +(dataframe["m"]*C_NaOH) 
                 #- V0*(BT/(1+dataframe["H"]/dataframe["KB"]))
                 - (V0)*(CO2/(1+(dataframe["H"]/(dataframe["K1_LK"]))+dataframe["K2_LK"]/dataframe["H"]))
                 - (V0)*(X1/(1+dataframe["H"]/K_X1))
                 #- (V0)*(X2/(1+dataframe["H"]/K_X3))
                 #- (V0)*(X2/(1+dataframe["H"]/K_X3))
                )


        return model - data_NaOH

    # create a set of Parameters
    params = Parameters()
    params.add('H0',   value= H0 ) #highest [H+] value used as initial estimate
    params.add('C_NaOH',   value= C_NaOH, vary = False ) 
    params.add('X1',   value= X1_initial, min = 0 )
    params.add('K_X1',   value= K_X1_initial, min = 0  )

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x_NaOH, data_NaOH))
    kws  = {'options': {'maxiter':100}}
    result = minner.minimize()

    # calculate final result
    final = data_NaOH + result.residual
    

    H0 = result.params.get('H0').value
    H0_std = result.params.get('H0').stderr

    X1 = result.params.get('X1').value
    X1_std = result.params.get('X1').stderr
    K_X1 = result.params.get('K_X1').value
    K_X1_std = result.params.get('K_X1').stderr
