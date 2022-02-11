#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 16:37:58 2022

@author: Dan Kerr & Charles Turner

This file will comprise the main class used in the sausage machine to perform 
organic alkalinity calculations.
"""

#########################
#IMPORT REQUIRED PACKAGES
#########################

import numpy as np
import pandas as pd
import math
import os
import matplotlib.pyplot as plt

from scipy.stats import linregress
from lmfit import Minimizer, Parameters, report_fit
from matplotlib.legend_handler import HandlerBase
from openpyxl import load_workbook
from IPython.display import Markdown, display

class org_alk_titration():
    def __init__(self,dataset_path,spreadsheet_name_TA = None
                                  ,spreadsheet_name_NaOH = None
                                  ,spreadsheet_name_BT = None):
        self.dataset_path = dataset_path
        self.spreadsheet_name_TA = spreadsheet_name_TA
        self.spreadsheet_name_NaOH = spreadsheet_name_NaOH
        self.spreadsheet_name_BT = spreadsheet_name_BT
        self.df_TA = None
        self.df_NaOH = None
        self.df_BT = None
        self.V0 = None 
        self.Va = None
        self.Vb = None
        self.V0_BT = None
        self.S_TA = None
        self.sample_id_TA = None
        self.data_start_TA = None
        self.initial_EV_TA = None
        self.initial_EV_BT = None
        self.initial_K_TA = None
        self.initial_K_BT = None
        self.initial_K_NaOh = None
        self.df_TA = None
        self.df_NaOH = None
        self.BT = None
        self.mass_TA_known = False
        self.mass_NaOH_known = False
        self.temp_TA_known = False
        self.TA_est_TA = None
        self.TA_est_BT = None
        self.E0_init_est_TA = None
        self.E0_init_est_BT = None
        self.H0 = None
        self.E0_final_TA = None
        self.E0_final_BT = None
        self.TA_final_TA = None
        self.TA_final_BT = None
        
    

    ## Magic Numbers ## 
    ###################
    #ACID TITRANT INFO# Needs to be user specified at some point hopefully
    ###################
    C = 0.10060392 # Estimate of OVERALL HCl titrant concentration (mol.kg-1) from CRM_182_0392 
    Cl_HCl = 0.10060392 #ionic strength of acid, 0.1M HCl made up in DI therefore [Cl-] = 0.1 M

    ########################################################
    #CONSTANTS CALCULATED BASED ON SALINITY AND TEMPERATURE# 
    ########################################################
    
    R = 8.314472 # Universal gas constant
    F = 96485.3399 # Faraday constant

    ###################
    #BASE TITRANT INFO# User specified too
    ###################
    C_NaOH = 0.082744091 # ± 0.000226775 determimed from acidic gran function of standardisation of NaOH using crm standardised HCl
    I_NaOH = 0.082744091 #ionic strength of NaOH 
    
    def read_excel_spreadsheets(self,TA_filename
                               ,NaOH_filename
                               ,BT_filename):
        # This function will read in the excel spreadsheets to memory 
        # containing the organic alkalinity titration
        
        TA_filename_full = os.path.join(self.dataset_path,TA_filename)
        NaOH_filename_full = os.path.join(self.dataset_path,NaOH_filename)
        BT_filename_full = os.path.join(self.dataset_path,BT_filename)
        
        self.df_TA = pd.read_excel(TA_filename_full)
        self.df_NaOH = pd.read_excel(NaOH_filename_full)
        self.df_BT = pd.read_excel(BT_filename_full)
        

    
    
    def extract_TA_data(self,start_idx=0,end_idx=9):
        # This should take the spreadsheet Dan gave me and save these data to 
        # the class instance. I've looked at it doesn't appear like the 
        # different titration classes can be all used as the same function.
        
        df_TA = self.df_TA
        
        self.V0 = df_TA.iloc[start_idx]['g_0']-df_TA.iloc[start_idx]['g_1'] # Sample mass (g)
        self.S_TA = df_TA.iloc[start_idx]['SALINITY']  # Sample Salinity 
        self.sample_id_TA = df_TA.iloc[start_idx]['SAMPLE']  # Sample ID
        self.data_start_TA = int(df_TA.iloc[start_idx]['data_start']-1) #row which titration starts, eg after initial acid addition and degassing
        self.initial_EV_TA = df_TA.iloc[end_idx]['102 Voltage (V)'] #EV of sample before any acid addition, at index = 10

    def extract_NaOH_data(self):
        kelvin_offset = 273.15
        self.df_NaOH["T"] = (self.df_NaOH["SAMPLE Temperature (°C)"]+kelvin_offset) #create column for sample temperature (KELVIN) at each titration point
        self.df_NaOH["NaOH_T"] = self.df_NaOH["NaOH Temperature (°C)"] #create colume for temperature (Degrees Celsius) of NaOH upon addition to cell 
        self.temp_TA_known = True 

    def extract_BT_data(self,start_idx=0):
        
        if self.mass_TA_known == False:
            raise AssertionError("Total Alkalinity mass must be known. Run convert_vol_to_mass on TA data first")
        if self.mass_NaOH_known == False:
            raise AssertionError("NaOH mass must be known. Run convert_vol_to_mass on NaOH data first")
        
        df_TA = self.df_TA
        df_NaOH = self.df_NaOH
        
        self.Va = df_TA['m'][df_TA.index[-1]] #DEFINE TOTAL MASS OF ACID ADDED DURING FIRST (TA) TITRATION
        self.Vb = df_NaOH['m'][df_NaOH.index[-1]] #DEFINE TOTAL MASS OF BASE ADDED DURING NAOH TITRATION
        self.V0_BT = (self.V0+self.Va+self.Vb) # Sample mass accounting for additions of acid and base (g) 
        self.sample_id_BT = self.df_BT.iloc[start_idx]['SAMPLE']  # Sample ID
        self.data_start_BT = int(self.df_BT.iloc[start_idx]['data_start']-1) #row which titration starts, eg after initial acid addition and degassing
        
    def strip_data(self,titration_label,start_idx=0,data_start=41):
        # Data start being 41 makes no sense and needs to be cleaned up into 
        # something (more?) sensible
        if titration_label == "TA":
            dataframe = self.df_TA
        elif titration_label == "NaOH":
            dataframe = self.df_NaOH
        elif titration_label == "BT":
            dataframe = self.df_BT
            data_start = self.data_start_BT
        else:
            raise ValueError("Dataframe label not recognised")
        
        dataframe['E(V)'] = dataframe.drop(dataframe.index[start_idx:data_start]
                                          ,axis=0)['102 Voltage (V)']
        dataframe['E(V)'] = dataframe['E(V)'].shift(-(data_start))
        dataframe['Sample_T'] = dataframe.drop(dataframe.index[start_idx:data_start]
                                              ,axis=0)['SAMPLE Temperature (°C)']
        dataframe['Sample_T'] = dataframe['Sample_T'].shift(-(data_start))
        dataframe['Acid_T'] = dataframe.drop(dataframe.index[start_idx:data_start]
                                            ,axis=0)['ACID Temperature (°C)']
        dataframe['Acid_T'] = dataframe['Acid_T'].shift(-(data_start))
        stripped_dataframe = dataframe[['E(V)', 'Sample_T', 'Acid_T', "mL"]].copy() #copy above variables to new "stripped_df"
        stripped_dataframe = stripped_dataframe.dropna() #remove NaN values from stripped_df
        
        if titration_label == "TA":
            self.df_TA = stripped_dataframe
        elif titration_label == "NaOH":
            self.df_NaOH = stripped_dataframe
        elif titration_label == "BT":
            self.df_BT = stripped_dataframe


    def density_curve_info(self,titration_label):
        if titration_label == "NaOH":
            slope_rho = -0.014702658
            intercept_rho = 1.27068345
        elif titration_label == "TA" or titration_label == "BT":
            slope_rho = -0.0008958
            intercept_rho = 1.02119193
        else:
            raise ValueError("Titration label not recognised")
        return slope_rho, intercept_rho
        
    def vol_to_mass(self,titration_label,initial_mL=0):
        
        if self.temp_TA_known == False:
            raise AssertionError("NaOH temperature not known. Run extract NaOH data first.")
        
        slope_rho, intercept_rho = self.density_curve_info(titration_label)
        
        if titration_label == "TA":
            dataframe = self.df_TA
            df_T_label = "Acid_T"   
            initial_mL = dataframe.iloc[0]["mL"] 
        elif titration_label == "NaOH":
            dataframe = self.df_NaOH
            df_T_label = "NaOH_T"
            initial_mL = 0
        elif titration_label == "BT":
            dataframe = self.df_BT
            df_T_label = "Acid_T"
            initial_mL = dataframe.iloc[0]["mLs
                                           "]
        else:
            raise ValueError("Dataframe label not recognised")

        dataframe["rho"] = (dataframe[df_T_label]*slope_rho+intercept_rho) # Density of titrant g.cm-3 

           
        initial_g = initial_mL*dataframe.iloc[0]["rho"]#convert inital volume of base/acid added to mass value (g)
        dataframe['delta_mL'] = dataframe['mL'].diff() #calculate the incremental values of volume of base/acid added
        dataframe['delta_g'] = dataframe['delta_mL']*dataframe["rho"] #convert the incremental values of volume of base/acid added to mass values
        dataframe = dataframe.fillna(0) #initial value of df['delta_g'] will be NA by default, replace with 0
        dataframe['delta_g'] = dataframe['delta_g'].fillna(0)
        dataframe['m'] = initial_g+np.cumsum(dataframe['delta_g'])#calculate cumulative total of mass of base/acid added, (initial mass of acid added + increment mass 1,  initial mass of acid added + increment mass 1 + increment mass 2 ...)
        
        if titration_label == "TA":
            self.df_TA = dataframe
            self.mass_TA_known = True
            self.Va = dataframe['m'][dataframe.index[-1]]
        elif titration_label == "NaOH":
            self.df_NaOH = dataframe
            self.Vb = dataframe['m'][dataframe.index[-1]]
            self.mass_NaOH_known = True
        elif titration_label == "BT":
            self.df_BT = dataframe
        
    def nernst_factor(self,titration_label):
        
        if titration_label == "TA":
            dataframe = self.df_TA
        elif titration_label == "NaOH":
            dataframe = self.df_NaOH
        elif titration_label == "BT":
            dataframe = self.df_BT
        else:
            raise ValueError("Dataframe label not recognised")
        
        dataframe["T"] = (dataframe["Sample_T"]+273.15)  # CREATE COLUMN SAMPLE TEMPERATURE (KELVIN) AT EACH TITRATION POINT
        dataframe["K"] = (self.R*dataframe['T'])/self.F # Nernst factor 
        initial_K = dataframe.iloc[9]['K'] # Initial Nernst factor, used to calculate initial pH

        if titration_label == "TA":
            self.df_TA = dataframe
            self.initial_K_TA = initial_K
        elif titration_label == "NaOH":
            self.df_NaOH = dataframe
            self.initial_K_NaOH = initial_K
        elif titration_label == "BT":
            self.df_BT = dataframe
            self.initial_K_BT = initial_K
            self.initial_EV_BT = dataframe.iloc[0]['E(V)'] #EV of sample before any acid addition

            
    def ion_strength_salinity(self,titration_label):
        
        if titration_label == "TA":
            dataframe = self.df_TA
            S = self.S_TA
            V0 = self.V0
            titration_soln = self.Cl_HCl 
        elif titration_label == "NaOH":
            dataframe = self.df_NaOH
            df_TA = self.df_TA
            S = df_TA['S'][df_TA.index[-1]]
            V0 = self.V0 + self.Va 
            titration_soln = self.I_NaOH 
        elif titration_label == "BT":
            dataframe = self.df_BT
            df_NaOH = self.df_NaOH
            S = df_NaOH['S'][df_NaOH.index[-1]] # This is giving NaN. I suspect all problems stem from ehre=
            V0 = self.V0_BT
            titration_soln = self.Cl_HCl 

        ImO = (19.924*S/(1000-1.005*S))

        dataframe['ImO'] = (ImO*V0 + dataframe['m']*titration_soln)/(V0 + dataframe['m']) 
        dataframe['S'] = (1000*dataframe['ImO']) / (19.924 + 1.005*dataframe['ImO'])
      
        if titration_label == "TA":
            self.df_TA = dataframe
        elif titration_label == "NaOH":
            self.df_NaOH = dataframe
        elif titration_label == "BT":
            self.df_BT = dataframe
            
    def equilibrium_consts_sulfate_HF(self,titration_label):
        # Needs to be done after calculating ionic strength and salinity (same 
        # for TA and BT (similar has to be done for NaOH titration, bells &
        # whistles))
        if titration_label == "TA":
            dataframe = self.df_TA
        elif titration_label == "BT":
            dataframe = self.df_BT
        elif titration_label == "NaOH":
            # Piss off this bridge when we get to it
            x = 0
        else:
            raise ValueError("Dataframe label not recognised")
        
        
        dataframe['K_S'] = np.exp(-4276.1/dataframe['T'] + 141.328 
                                  - 23.093*np.log(dataframe['T'])
                                  +(-13856/dataframe['T'] + 324.57- 47.986*np.log(dataframe['T']))*
                                  (dataframe['ImO']**(1/2)) 
                                  +(35474/dataframe['T']  - 771.54 +114.723*np.log(dataframe['T']))
                                  *dataframe['ImO'] - (2698/dataframe['T'])*(dataframe['ImO'])**(3/2)+(1776/dataframe['T'])
                                  *(dataframe['ImO'])**(2) +np.log(1-0.001005*dataframe['S'])) # pKa Bisulfate ion [HSO4-] K_S from Dickson1990
        
        
        dataframe['K_F'] = np.exp(874/dataframe['T'] - 9.68 + 0.111*dataframe['S']**0.5)  # pKa Hydrogen Fluoride ion [HF] K_F  from Perex & Fraga 1987
        
        dataframe['S_T'] =  (0.14/96.062)*(dataframe['S']/1.80655)# Total Sulphate concentration S_T from Morris & Riley 1966
        dataframe['F_T'] = (0.000067/18.998)*(dataframe['S']/1.80655) # Total Hydrogen Fluoride concentration F_T from Riley 1996
        dataframe["Z"] = (1+dataframe['S_T']/dataframe['K_S']) #from dickson 2007 Z = (1+S_T/K_S)

        if titration_label == "TA":
            self.df_TA = dataframe
        elif titration_label == "BT":
            self.df_BT = dataframe
            
    def gran_func(self,titration_label):
        
        if titration_label == "TA":
            dataframe = self.df_TA
            V0 = self.V0
        elif titration_label == "BT":
            dataframe = self.df_BT
            V0 = self.V0_BT
        else:
            raise ValueError("Dataframe label not recognised")
        
        dataframe["F1"] = ((V0+dataframe["m"])*np.exp((dataframe["E(V)"]/(dataframe['K'])))) #Calculate Gran Function F1 at each titration point
        dataframe = dataframe[dataframe["F1"].between(10000, 1000000)] #drop all gran funciton values less than 10000 as these are TYPICALLY non linear with respect to m
        slope, intercept, r_value, p_value, std_err =linregress(dataframe["m"], dataframe["F1"])#CALL SLOPE AND INTERCEPT OF Gran function F1
        equivalence_point = -intercept/slope #Calculate equivalence point estimate (g) from Gran function F1
        TA_est = (equivalence_point*self.C)/V0 #Estimate TA using equivalence point estimate
        
        dataframe["E0_est"] = (dataframe["E(V)"]-(dataframe["K"])*np.log((-V0*TA_est + dataframe["m"]*self.C)/(V0 + dataframe["m"]))) #CALCULATE EO ESTIMATE FOR EACH TITRATION POINT
        E0_init_est = dataframe["E0_est"].mean()#AVERAGE ALL EO VALUES TO OBTAIN AN INITIAL ESTIMATE OF EO
        dataframe["H"] = (np.exp((dataframe["E(V)"]-E0_init_est)/(dataframe["K"]))) #USING EO INITIAL ESTIMATE CALCULATE [H'] FOR EACH TITRATION POINT
        dataframe["GRAN_pH"] = -(np.log10(np.exp((dataframe["E(V)"]-E0_init_est)/dataframe["K"])))#CALCULATE GRAN pH

        if titration_label == "TA":
            self.df_TA = dataframe
            self.TA_est_TA = TA_est
            self.E0_init_est_TA = E0_init_est
        elif titration_label == "BT":
            self.df_BT = dataframe
            self.TA_est_BT = TA_est
            self.E0_init_est_BT = E0_init_est
            
    def nl_least_squares(self,titration_label):
        
        if titration_label == "TA":
            dataframe = self.df_TA
            E0_init_est = self.E0_init_est_TA
            initial_EV = self.initial_EV_TA
            TA_est = self.TA_est_TA
            initial_K = self.initial_K_TA
            V0 = self.V0
        elif titration_label == "BT":
            dataframe = self.df_BT
            E0_init_est = self.E0_init_est_BT
            initial_EV = self.initial_EV_BT
            TA_est = self.TA_est_BT
            initial_K = self.initial_K_BT
            V0 = self.V0_BT
        else:
            raise ValueError("Dataframe label not recognised")
            
        new_dataframe = dataframe[dataframe["GRAN_pH"].between(3, 3.5)]#SELECT ONLY TITRATION POINTS WHICH pH ARE BETWEEN 3.0 - 3.5 
        data_points = len(new_dataframe.index) #CALCULATE NUMBER OF DATA POINTS IN 3-3.5 pH RANGE

        # DEFINE FUNCTION WHICH RETURNS VALUE(S) TO BE MINIMISED, IN THIS CASE SSR 
        x = new_dataframe["H"]
        data = new_dataframe["m"]

        def fcn2min(params, x, data):
            f_NLSF = params['f_NLSF']
            TA_est_NLSF = params['TA_est_NLSF']
            model = ((np.sum((TA_est_NLSF + 
                              ((V0+new_dataframe["m"])/V0)*
                              ((f_NLSF*new_dataframe["H"])/new_dataframe["Z"]) -
                              (new_dataframe["m"]/V0)*self.C)**2))*10**12)
    
            return model - data
    
        # create a set of Parameters
        params = Parameters()
        params.add('f_NLSF',   value= 1)
        params.add('TA_est_NLSF', value= TA_est)
        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(x, data))
        kws  = {'options': {'maxiter':10}}
        result = minner.minimize()
        # calculate final result
        final = data + result.residual                                                                        
        
        
        #################################
        #EXTRACT AND PROCESS TA NLSF RESULTS
        #################################    
        TA_processed = result.params.get('TA_est_NLSF').value #EXTRACT INTIAL TA VALUE M/KG-1
        TA_final = TA_processed*10**6 #COVERT INTIAL TA VALUE TO µmol/KG-1
        f = result.params.get('f_NLSF').value #EXTRACT NLSF F VALUE                               
        E0_processed = E0_init_est + dataframe["K"]*np.log(f) #CALCULATE E0 FROM NLSF F VALUE               
        E0_final = E0_processed.mean() #FINAL ESTIMATE OF E0
        new_dataframe["pH"] = -np.log10(np.exp((dataframe["E(V)"]-E0_final)/dataframe["K"])) #CALCULATE pH AT EACH TITRATION POINT FROM E0 FINAL ESTIMATE                                                                          
        initial_pH = -np.log10(np.exp((initial_EV-E0_final)/initial_K)) #CALCULATE INTIAL pH FROM E0 FINAL ESTIMATE
        
        if titration_label == "TA":
            self.E0_final_TA = E0_final
            self.TA_final_TA = TA_final
        elif titration_label == "BT":
            self.E0_final_BT = E0_final
            self.TA_final_BT = TA_final
            
        print(titration_label,": ", TA_final, 'µmol.kg-1')
        
    def pH_H_OH_H0_conc(self):
        dataframe = self.df_NaOH        
        
        dataframe["K"] = (self.R*dataframe["T"])/self.F #Get K value at eachpoint during NaOH titration
        dataframe["pH"] = -np.log10(np.exp((dataframe["102 Voltage (V)"]-self.E0_final_TA)/dataframe["K"]))#Using EO estimated from TA NLSF procedure calculate pH at each point during NaOH titration
        dataframe["H"] = 10**-(dataframe["pH"]) #Using pH calculate [H+] at each titration point
        dataframe['pKw'] = (-np.log10(np.exp(148.9652-13847.26/dataframe["T"] -
                           23.6521*np.log(dataframe["T"])+(-5.977+118.67/dataframe["T"] + 
                           1.0495*np.log(dataframe["T"]))*dataframe["S"]**0.5-0.01615*dataframe["S"]))) #Acid dissociation constant of Water
        dataframe['OH'] = (10**-dataframe["pKw"])/dataframe["H"] #using Acid dissociation constant of Water and [H+] calculate [OH-]
        initial_EV_NaOH = dataframe.iloc[0]['102 Voltage (V)'] #Proton concentration prior to NaOH addition, H0, defined as [H+] at end of first (TA) titration (start of NaOH titration)
        initial_K_NaOH = dataframe.iloc[0]['K']
        self.H0 = (np.exp((initial_EV_NaOH-self.E0_final_TA)/(initial_K_NaOH)))
        
        
    def pipeline(self):
        self.extract_TA_data()
        self.strip_data("TA")
        self.nernst_factor("TA")
        self.extract_NaOH_data()
        self.vol_to_mass("TA")
        self.vol_to_mass("NaOH")
        self.ion_strength_salinity("TA")
        self.equilibrium_consts_sulfate_HF("TA")
        self.gran_func("TA")
        self.nl_least_squares("TA")


        self.ion_strength_salinity("NaOH")
        self.pH_H_OH_H0_conc()

        self.extract_BT_data()
        self.strip_data("BT")
        self.vol_to_mass("BT")
        self.nernst_factor("BT")
        self.ion_strength_salinity("BT")
        self.equilibrium_consts_sulfate_HF("BT")
        self.gran_func("BT")
        self.nl_least_squares("BT")






