#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:47:18 2022

@author: ct
"""

import SausageMachine


tit = SausageMachine.org_alk_titration("~/Python/org-alk-sausage-machine")

tit.read_excel_spreadsheets("01.09.21.50UM.001_PROCESSED.xlsx",
                                       "01.09.21.50UM.001.NAOH_PROCESSED.xlsx",
                                       "01.09.21.50UM.001.BT_PROCESSED.xlsx")

tit.pipeline()

tit.init_minimiser()

tit.dissociation_consts()

# This is where it goes tits up 
tit.minimise(1)

