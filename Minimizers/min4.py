##################################################################
############################## 004 ###############################
##################################################################
# H0 fixed
# X1, KX1 fixed
# All data used to constrain X2, KX2  % Is this description correct?
##################################################################

def minner_004(dataframe):
    global H0, X1
    global C_NaOH, K_X1, X2, K_X2
    global X3, K_X3

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
        X3 = params['X3']
        K_X3 = params['K_X3']

        model = ((V0 + Va+ dataframe["m"])*(dataframe["H"]-dataframe["OH"]) 
                 -((V0+Va)*H0)
                 +(dataframe["m"]*C_NaOH) 
                  #- V0*(BT/(1+dataframe["H"]/dataframe["KB"]))
                 - (V0)*(CO2/(1+(dataframe["H"]/(dataframe["K1_LK"]))+dataframe["K2_LK"]/dataframe["H"]))
                 - (V0)*(X1/(1+dataframe["H"]/K_X1))
                 - (V0)*(X2/(1+dataframe["H"]/K_X2))
                 - (V0)*(X3/(1+dataframe["H"]/K_X3))
                )


        return model - data_NaOH

    # create a set of Parameters
    params = Parameters()
    params.add('H0',   value= H0, vary=False ) #highest [H+] value used as initial estimate
    params.add('C_NaOH',   value= C_NaOH, vary=False ) 
    params.add('X1',   value= X1 , vary=False )
    params.add('K_X1',   value= K_X1, vary=False )
    params.add('X2',   value= X2, vary=False )
    params.add('K_X2',   value= K_X2, vary=False )
    params.add('X3',   value= X3)
    params.add('K_X3',   value= K_X3)

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x_NaOH, data_NaOH))
    kws  = {'options': {'maxiter':100}}
    result = minner.minimize()
    
    # calculate final result
    final = data_NaOH + result.residual

    X3 = result.params.get('X3').value
    X3_std = result.params.get('X3').stderr

    K_X3 = result.params.get('K_X3').value
    K_X3_std = result.params.get('K_X3').stderr
