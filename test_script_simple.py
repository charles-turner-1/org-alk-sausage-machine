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

tit.minimise(1)

import SausageMachine

tit_simple = SausageMachine.org_alk_titration("~/Python/org-alk-sausage-machine")

tit_simple.read_excel_spreadsheets("01.09.21.50UM.001_PROCESSED.xlsx",
                                       "01.09.21.50UM.001.NAOH_PROCESSED.xlsx",
                                       "01.09.21.50UM.001.BT_PROCESSED.xlsx")

tit_simple.pipeline()
tit_simple.init_minimiser()
tit_simple.dissociation_consts()

k = tit_simple.minimise_01()

# So where the fuck is the error coming from?