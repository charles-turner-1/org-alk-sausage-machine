##################################################################
############################## 001 ###############################
##################################################################
# Data in pH range 0-5 used to constrain H0, X1, K_X1
##################################################################

def minner_001():
    global H0, X1, K_X1

    x_NaOH = df_constrain_001.m
    data_NaOH = df_constrain_001.H


    # DEFINE FUNCTION WHICH RETURNS VALUE(S) TO BE MINIMISED, IN THIS CASE SSR 
    def fcn2min(params, x_NaOH, data_NaOH):
        H0 = params['H0']
        C_NaOH = params['C_NaOH']
        X1 = params['X1']
        K_X1 = params['K_X1']

        model = ((V0 + Va+ df_constrain_001["m"])*(df_constrain_001["H"]-df_constrain_001["OH"]) 
                 -((V0+Va)*H0)
                 +(df_constrain_001["m"]*C_NaOH) 
                 #- V0*(BT/(1+df_constrain_001["H"]/df_constrain_001["KB"]))
                 - (V0)*(CO2/(1+(df_constrain_001["H"]/(df_constrain_001["K1_LK"]))+df_constrain_001["K2_LK"]/df_constrain_001["H"]))
                 - (V0)*(X1/(1+df_constrain_001["H"]/K_X1))
                 #- (V0)*(X2/(1+df_constrain_001["H"]/K_X3))
                 #- (V0)*(X2/(1+df_constrain_001["H"]/K_X3))
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
