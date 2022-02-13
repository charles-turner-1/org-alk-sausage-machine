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



##################################################################
############################## 002 ###############################
##################################################################
# H0 fixed
# Data in pH range 0-7.5 to fix X1, KX1
##################################################################

def minner_002():
    global H0, C_NaOH, X1, K_X1, X2, K_X2

    x_NaOH = df_constrain_002.m
    data_NaOH = df_constrain_002.H
    # DEFINE FUNCTION WHICH RETURNS VALUE(S) TO BE MINIMISED, IN THIS CASE SSR 
    def fcn2min(params, x_NaOH, data_NaOH):
        H0 = params['H0']
        C_NaOH = params['C_NaOH']
        X1 = params['X1']
        K_X1 = params['K_X1']
        X2 = params['X2']
        K_X2 = params['K_X2']

        model = ((V0 + Va+ df_constrain_002["m"])*(df_constrain_002["H"]-df_constrain_002["OH"]) 
                 -((V0+Va)*H0)
                 +(df_constrain_002["m"]*C_NaOH) 
                 #- V0*(BT/(1+df_constrain_002["H"]/df_constrain_002["KB"]))
                 - (V0)*(CO2/(1+(df_constrain_002["H"]/(df_constrain_002["K1_LK"]))+df_constrain_002["K2_LK"]/df_constrain_002["H"]))
                 - (V0)*(X1/(1+df_constrain_002["H"]/K_X1))
                 - (V0)*(X2/(1+df_constrain_002["H"]/K_X2))
                 #- (V0)*(X2/(1+df_constrain_002["H"]/K_X3))
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



##################################################################
############################## 003 ###############################
##################################################################
# H0 fixed
# X1, KX1 fixed
# All data used to constrain X2, KX2
##################################################################

def minner_003():
    global H0, C_NaOH, X1, K_X1, X2, K_X2, X3, K_X3
    
    x_NaOH = df_NaOH.m
    data_NaOH = df_NaOH.H
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

        model = ((V0 + Va+ df_NaOH["m"])*(df_NaOH["H"]-df_NaOH["OH"]) 
                 -((V0+Va)*H0)
                 +(df_NaOH["m"]*C_NaOH) 
                 #- V0*(BT/(1+df_NaOH["H"]/df_NaOH["KB"]))
                 - (V0)*(CO2/(1+(df_NaOH["H"]/(df_NaOH["K1_LK"]))+df_NaOH["K2_LK"]/df_NaOH["H"]))
                 - (V0)*(X1/(1+df_NaOH["H"]/K_X1))
                 - (V0)*(X2/(1+df_NaOH["H"]/K_X2))
                 - (V0)*(X3/(1+df_NaOH["H"]/K_X3))
                )


        return model - data_NaOH

    # create a set of Parameters
    params = Parameters()
    params.add('H0',   value= H0, vary=False ) #highest [H+] value used as initial estimate
    params.add('C_NaOH',   value= C_NaOH, vary=False ) 
    params.add('X1',   value= X1 , vary=False )
    params.add('K_X1',   value= K_X1, vary=False )
    params.add('X2',   value= X2)
    params.add('K_X2',   value= K_X2 )
    params.add('X3',   value= X3_initial, min = 0 )
    params.add('K_X3',   value= K_X3_initial, min = 0 )

    # do fit, here with leastsq model
    minner = Minimizer(fcn2min, params, fcn_args=(x_NaOH, data_NaOH))
    kws  = {'options': {'maxiter':100}}
    result = minner.minimize()
    
    # calculate final result
    final = data_NaOH + result.residual

    H0 = result.params.get('H0').value
    H0_std = result.params.get('H0').stderr


    C_NaOH = result.params.get('C_NaOH').value 

    X2 = result.params.get('X2').value
    X2_std = result.params.get('X2').stderr

    K_X2 = result.params.get('K_X2').value
    K_X2_std = result.params.get('K_X2').stderr

    X3 = result.params.get('X3').value
    X3_std = result.params.get('X3').stderr

    K_X3 = result.params.get('K_X3').value
    K_X3_std = result.params.get('K_X3').stderr
    
##################################################################
############################## 004 ###############################
##################################################################
# H0 fixed
# X1, KX1 fixed
# All data used to constrain X2, KX2
##################################################################

def minner_004():
    global H0, C_NaOH, X1, K_X1, X2, K_X2, X3, K_X3

    x_NaOH = df_NaOH.m
    data_NaOH = df_NaOH.H
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

        model = ((V0 + Va+ df_NaOH["m"])*(df_NaOH["H"]-df_NaOH["OH"]) 
                 -((V0+Va)*H0)
                 +(df_NaOH["m"]*C_NaOH) 
                  #- V0*(BT/(1+df_NaOH["H"]/df_NaOH["KB"]))
                 - (V0)*(CO2/(1+(df_NaOH["H"]/(df_NaOH["K1_LK"]))+df_NaOH["K2_LK"]/df_NaOH["H"]))
                 - (V0)*(X1/(1+df_NaOH["H"]/K_X1))
                 - (V0)*(X2/(1+df_NaOH["H"]/K_X2))
                 - (V0)*(X3/(1+df_NaOH["H"]/K_X3))
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
