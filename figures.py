def figure_001():


    df_constrain_001["m_calc_001"] = ( (CO2*(V0)/((df_constrain_001["H"]/df_constrain_001["K1_LK"])+(df_constrain_001["K2_LK"]/df_constrain_001["H"])+1) 
                                         + X1*(V0)/((df_constrain_001["H"]/K_X1)+1)
                                        #+ X2*(V0)/((df_constrain_001["H"]/K_X2)+1)
                                        #+ X3*(V0)/((df_constrain_001["H"]/K_X3)+1)
                                        #+ V0*(BT/(1+df_constrain_001["H"]/df_constrain_001["KB"]))
                                         + H0*(V0+Va)
                                         -V0*(df_constrain_001["H"]-df_constrain_001["OH"])
                                         -Va*(df_constrain_001["H"]-df_constrain_001["OH"]))/(df_constrain_001["H"]-df_constrain_001["OH"]+C_NaOH) )

    ssr_001 = np.sum((df_constrain_001['m']-df_constrain_001['m_calc_001'])**2)
    ssr_001_trunc = truncate(ssr_001, 5)
    
    #KO2016 model CO2 term
    x_meas = df_constrain_001["m"]
    x_calc = df_constrain_001["m_calc_001"]
    plt.xlabel('NaOH added (g)', fontsize=18)
    y_meas = df_constrain_001["pH"]
    y_calc = df_constrain_001["pH"]
    plt.ylabel('pH', fontsize=16)
    #plt.axhline(0, color='red', linestyle='--') 
    #plt.axhline(0, color='red', linestyle='--') 
    graph = plt.scatter(x_meas, y_meas, c = 'black', marker = "1")
    graph = plt.plot(x_calc, y_calc, c = 'red')
    plt.grid(False)
    ax = plt.gca()
    ax.tick_params(bottom='on', left='on', labelleft='on', labelbottom='on', length=5, labelsize = 10.5)
    plt.rc('axes',edgecolor='black')
    plt.annotate('SSR: {}'.format(ssr_001_trunc), xy=(0.0650, 0.75), xycoords='axes fraction')


    list_color  = ["black","red",]
    list_mak    = ["1",       "_"]
    list_lab    = ['Measured','Calculated',]

    ax.legend(list(zip(list_color,list_mak)), list_lab, 
              handler_map={tuple:MarkerHandler()}) 

    plt.show()
    
    print('X1 (initial):', X1*10**6, "| pK1(initial): ", -np.log10(K_X1), '| H0 :', H0 ) 

def figure_002(): 
    global ssr_002
    
        #002

    df_constrain_002["m_calc_002"] = ( (CO2*(V0)/((df_constrain_002["H"]/df_constrain_002["K1_LK"])+(df_constrain_002["K2_LK"]/df_constrain_002["H"])+1) 
                                         + X1*(V0)/((df_constrain_002["H"]/K_X1)+1)
                                        + X2*(V0)/((df_constrain_002["H"]/K_X2)+1)
                                        #+ X3*(V0)/((df_constrain_002["H"]/K_X3)+1)
                                        #+ V0*(BT/(1+df_constrain_002["H"]/df_constrain_002["KB"]))
                                         + H0*(V0+Va)
                                         -V0*(df_constrain_002["H"]-df_constrain_002["OH"])
                                         -Va*(df_constrain_002["H"]-df_constrain_002["OH"]))/(df_constrain_002["H"]-df_constrain_002["OH"]+C_NaOH) )

    ssr_002 = np.sum((df_constrain_002['m']-df_constrain_002['m_calc_002'])**2)
    ssr_002_trunc = truncate(ssr_002, 5)


    
    #KO2016 model CO2 term
    x_meas = df_constrain_002["m"]
    x_calc = df_constrain_002["m_calc_002"]
    plt.xlabel('NaOH added (g)', fontsize=18)
    y_meas = df_constrain_002["pH"]
    y_calc = df_constrain_002["pH"]
    plt.ylabel('pH', fontsize=16)
    graph = plt.scatter(x_meas, y_meas, c = 'black', marker = "1")
    graph = plt.plot(x_calc, y_calc, c = 'red')
    plt.grid(False)
    ax = plt.gca()
    ax.tick_params(bottom='on', left='on', labelleft='on', labelbottom='on', length=5, labelsize = 10.5)
    plt.rc('axes',edgecolor='black')
    plt.annotate('SSR: {}'.format(ssr_002_trunc), xy=(0.0650, 0.75), xycoords='axes fraction')


    list_color  = ["black","red",]
    list_mak    = ["1",       "_"]
    list_lab    = ['Measured','Calculated',]


    plt.show()


    print('X1:', X1*10**6, "| pK1: ", -np.log10(K_X1), '| Deviation % (g) NaOH :', (ssr_002/Vb)*100 ) 
    
    
def figure_003():  
    global ssr_003
    
        #003

    df_NaOH["m_calc_003"] = ( (CO2*(V0)/((df_NaOH["H"]/df_NaOH["K1_LK"])+(df_NaOH["K2_LK"]/df_NaOH["H"])+1) 
                                         + X1*(V0)/((df_NaOH["H"]/K_X1)+1)
                                        + X2*(V0)/((df_NaOH["H"]/K_X2)+1)
                                        + X3*(V0)/((df_NaOH["H"]/K_X3)+1)
                                         #+V0*(BT/(1+df_NaOH["H"]/df_NaOH["KB"]))
                                         + H0*(V0+Va)
                                         -V0*(df_NaOH["H"]-df_NaOH["OH"])
                                         -Va*(df_NaOH["H"]-df_NaOH["OH"]))/(df_NaOH["H"]-df_NaOH["OH"]+C_NaOH) )

    ssr_003 = np.sum((df_NaOH['m']-df_NaOH['m_calc_003'])**2)
    ssr_003_trunc = truncate(ssr_003, 5)
    
    
    #KO2016 model CO2 term
    x_meas = df_NaOH["m"]
    x_calc = df_NaOH["m_calc_003"]
    plt.xlabel('NaOH added (g)', fontsize=18)
    y_meas = df_NaOH["pH"]
    y_calc = df_NaOH["pH"]
    plt.ylabel('pH', fontsize=16)
    graph = plt.scatter(x_meas, y_meas, c = 'black', marker = "1")
    graph = plt.plot(x_calc, y_calc, c = 'red')
    plt.grid(False)
    ax = plt.gca()
    ax.tick_params(bottom='on', left='on', labelleft='on', labelbottom='on', length=5, labelsize = 10.5)
    plt.rc('axes',edgecolor='black')
    plt.annotate('SSR: {}'.format(ssr_003_trunc), xy=(0.0650, 0.75), xycoords='axes fraction')


    list_color  = ["black","red",]
    list_mak    = ["1",       "_"]
    list_lab    = ['Measured','Calculated',]

    plt.show()


    print('X2:', X2*10**6, "| pK2: ", -np.log10(K_X2), '| Deviation % (g) NaOH :', (ssr_003/Vb)*100 ) 
    
    
def figure_004(): 
    global ssr_004
    
        #004

    df_NaOH["m_calc_004"] = ( (CO2*(V0)/((df_NaOH["H"]/df_NaOH["K1_LK"])+(df_NaOH["K2_LK"]/df_NaOH["H"])+1) 
                                         + X1*(V0)/((df_NaOH["H"]/K_X1)+1)
                                        + X2*(V0)/((df_NaOH["H"]/K_X2)+1)
                                        + X3*(V0)/((df_NaOH["H"]/K_X3)+1)
                                         + H0*(V0+Va)
                                         #+V0*(BT/(1+df_NaOH["H"]/df_NaOH["KB"]))
                                         -V0*(df_NaOH["H"]-df_NaOH["OH"])
                                         -Va*(df_NaOH["H"]-df_NaOH["OH"]))/(df_NaOH["H"]-df_NaOH["OH"]+C_NaOH) )

    ssr_004 = np.sum((df_NaOH['m']-df_NaOH['m_calc_004'])**2)
    ssr_004_trunc = truncate(ssr_004, 5)
        #KO2016 model CO2 term
    x_meas = df_NaOH["m"]
    x_calc = df_NaOH["m_calc_004"]
    plt.xlabel('NaOH added (g)', fontsize=18)
    y_meas = df_NaOH["pH"]
    y_calc = df_NaOH["pH"]
    plt.ylabel('pH', fontsize=16)
    graph = plt.scatter(x_meas, y_meas, c = 'black', marker = "1")
    graph = plt.plot(x_calc, y_calc, c = 'red')
    plt.grid(False)
    ax = plt.gca()
    ax.tick_params(bottom='on', left='on', labelleft='on', labelbottom='on', length=5, labelsize = 10.5)
    plt.rc('axes',edgecolor='black')
    plt.annotate('SSR: {}'.format(ssr_004_trunc), xy=(0.0650, 0.75), xycoords='axes fraction')


    list_color  = ["black","red",]
    list_mak    = ["1",       "_"]
    list_lab    = ['Measured','Calculated',]


    ax.legend(list(zip(list_color,list_mak)), list_lab, 
              handler_map={tuple:MarkerHandler()}) 

    plt.show()


    print('X3:', X3*10**6, "| pK3: ", -np.log10(K_X3), '| Deviation % (g) NaOH :', (ssr_004/Vb)*100 ) 
