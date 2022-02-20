#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:47:18 2022

@author: ct
"""

import SausageMachine


titration_test = SausageMachine.org_alk_titration("/home/ct/DanKerr/")

titration_test.read_excel_spreadsheets("01.09.21.50UM.001_PROCESSED.xlsx",
                                       "01.09.21.50UM.001.NAOH_PROCESSED.xlsx",
                                       "01.09.21.50UM.001.BT_PROCESSED.xlsx")

titration_test.pipeline()




# tit = SausageMachine.org_alk_titration("/home/ct/DanKerr")
# tit.read_excel_spreadsheets("01.09.21.50UM.001_PROCESSED.xlsx",
#                             "01.09.21.50UM.001.NAOH_PROCESSED.xlsx",
#                             "01.09.21.50UM.001.BT_PROCESSED.xlsx")


# tit.extract_TA_data() # Unsuprisingly, no errors here
# tit.strip_data("TA")
# tit.nernst_factor("TA")
# tit.extract_NaOH_data()
# tit.vol_to_mass("TA")