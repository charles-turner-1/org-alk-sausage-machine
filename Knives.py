#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 16:40:59 2022

@author: Dan Kerr & Charles Turner

This file contains additional functions contained within the organic alkalinity
sausage machine that I (CT) have decreed are tangential enough that they should 
be contained in a separate file.
"""

def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper

def printmd(string):
    display(Markdown(string))
    
class MarkerHandler(HandlerBase):
    def create_artists(self, legend, tup,xdescent, ydescent,
                        width, height, fontsize,trans):
        return [plt.Line2D([width/2], [height/2.],ls="",
                       marker=tup[1],color=tup[0], transform=trans)]
