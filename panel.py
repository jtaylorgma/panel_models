########################################################################################################################
#Program: panel.py
#Project: Panel Models
#Author: Josh Taylor
#Last Edited: 6/18/15
########################################################################################################################

""" Estimate panel data models

    created: 6/18/15


    Notes
    -----
   This is only intended to be a very rudimentary module for implementing panel data models. There will be many things
   that this module cannot do, but hopefully this model will be sufficient to accomplish basic tasks using panel models
   until someone else can come up with something better.
"""

import statsmodels.api as sm
import pandas as pd
import numpy as np


class PanelModel:
    def __init__(self, formula, effects = "random", robust = False, data):
        legalEffects = ["fixed", "random"]
        self.formula = formula
        if effects in legalEffects:
            self.effects = effects
        else:
            print "Effects must be: "
            for i in range(len(legalEffects)):
                print
                print legalEffects[i]
                print
                if i < len(legalEffects) - 1:
                    print "or"
            raise ValueError("could not implement panel model")
        if robust != False and robust != True:
            print
            raise ValueError("robust parameter can only take a boolean value")
        self.robust = robust
        if type(data) is not pd.core.panel.Panel:
            raise ValueError("Data must be a panel dataframe, i.e. of type pandas.core.panel.Panel")


    def balanceChecker(self):


    def fit(self):
        if self.effects == "fixed":
            pass
        elif self.effects == 'random':
            self.randomEffects()

    def randomEffects(self):
