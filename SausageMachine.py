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
        self.S_TA = None
        self.sample_id_TA = None
        self.data_start_TA = None
        self.initial_EV_TA = None
        self.df_TA = None
        self.df_NaOH = None
        self.BT = None
        
    

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
    
    def extract_BT_data(start_idx=0):
        self.V0_BT = (self.V0+Va+Vb) # Sample mass accounting for additions of acid and base (g) 
        self.sample_id_BT = df_BT.iloc[start_idx]['SAMPLE']  # Sample ID
        self.data_start_BT = int(df_BT.iloc[start_idx]['data_start']-1) #row which titration starts, eg after initial acid addition and degassing
        
    def strip_data(self,dataframe_name,start_idx=0,data_start=42):
        if dataframe_name == "TA":
            dataframe = self.df_TA
        elif dataframe_name == "NaOH":
            dataframe = self.df_NaOH
        elif dataframe_name == "BT":
            dataframe = self.df_BT
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
        
        if dataframe_name == "TA":
            self.df_TA = stripped_dataframe
        elif dataframe_name == "NaOH":
            self.df_NaOH = stripped_dataframe
        elif dataframe_name == "BT":
            self.df_BT = stripped_dataframe


    def density_curve_info(titration_label):
        if titration_label == "NaOH":
            slope_rho = -0.014702658
            intercept_rho = 1.27068345
        elif titration_label == "TA":
            slope_rho = -0.0008958
            intercept_rho = 1.02119193
        else:
            raise ValueError("Titration label not recognised")
        return slope_rho, intercept_rho
        
    def convert_vol_to_mass(titration_label,dataframe,initial_mL_base=0):
        slope_rho, intercept_rho = density_curve_info(titration_label)
        df_T_label = titration_label + "_T"

        dataframe["rho"] = (dataframe[df_T_label]*slope_rho+intercept_rho) # Density of titrant g.cm-3 
        initial_g_base = initial_mL_base*dataframe.iloc[0]["rho"]#convert inital volume of base added to mass value (g)
        dataframe['delta_mL'] = dataframe['mL'].diff() #calculate the incremental values of volume of base added
        dataframe['delta_g'] = dataframe['delta_mL']*dataframe["rho"] #convert the incremental values of volume of base added to mass values
        #dataframe = dataframe.fillna(0) #initial value of df['delta_g'] will be NA by default, replace with 0
        dataframe['delta_g'] = dataframe['delta_g'].fillna(0)
        dataframe['m'] = initial_g_base+np.cumsum(dataframe['delta_g'])#calculate cumulative total of mass of base added, (initial mass of acid added + increment mass 1,  initial mass of acid added + increment mass 1 + increment mass 2 ...)
        
    def calc_nernst_factor(titration_label,dataframe):
        T = dataframe["T"] = (dataframe["Sample_T"]+273.15)  # CREATE COLUMN SAMPLE TEMPERATURE (KELVIN) AT EACH TITRATION POINT
        dataframe["K"] = (R*T)/F # Nernst factor 
        initial_K = dataframe.iloc[9]['K'] # Initial Nernst factor, used to calculate initial pH
        if titration_label == "BT":
            initial_EV = dataframe.iloc[0]['E(V)'] #EV of sample before any acid addition

        




















