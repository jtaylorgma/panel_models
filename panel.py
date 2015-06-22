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

   Reference: Wooldridge, Econometric Analysis of Cross-Section and Panel Data
"""

import statsmodels.formula.api as sm
import pandas as pd
import numpy as np
import re


class PanelModel:
    def __init__(self, formula, effects = "random", time_fe = False, entity_fe = False, robust = False, data = None):
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
            raise ValueError("could not create panel model")
        if robust != False and robust != True:
            print
            raise ValueError("robust parameter can only take a boolean value")
        self.robust = robust
        if type(data) is not pd.core.panel.Panel:
            raise ValueError("Data must be a panel dataframe, i.e. of type pandas.core.panel.Panel")
        self.panel = data
        self.balanced = True
        self.timeVar = self.panel.minor_axis
        self.idVar = self.panel.major_axis
        self.indepVars = re.findall(r"\w+",formula)[1:]
        self.depVar = re.findall(r"\w+",formula)[0]

        if time_fe:
            self.formula = self.formula + "C(%s)" %self.timeVar
        if entity_fe and effects == "random":
            self.formula = self.formula + "C(%s)" %self.idVar

        ### attributes to be defined later:
        self.bse = None
        self.cov_params = None
        self.params = None
        self.feParams = None
        self.pvalues = None
        self.betas = None
        self.Sigma = None
        self.SigmaRobust = None
        self.omega = None


    def balanceChecker(self):
        message = "This panel is balanced"
        for id in self.idVar:
            for t in self.timeVar:
                if t not in self.panel.major_xs(self.panel.major_axis[id]).index: #checking if each id has all times
                    message = "This panel is unbalanced"
                    self.balanced = False
                    break
        print message

    def fit(self):
        if self.effects == "fixed":
            self.fixedEffects()
        elif self.effects == 'random':
            self.randomEffects()

    def randomEffects(self):
        self.balanceChecker()
        if self.balanced:
            self.betasREBalanced()
        else:
            if self.robust:
                pass
            else:
                pass

    def fixedEffects(self):
        pooledDf = self.panel.to_frame().reset_index()
        formulaFE = self.formula + " + C(%s)" %self.idVar
        if not self.robust:
            print "Sorry for the inconvenience but this result is of type: statsmodels.regression.linear_model.RegressionResultsWrapper"
        FE = sm.ols(formula=formulaFE, data = pooledDf).fit()
        if self.robust:
            pooledDf['resid'] = FE.resid
            groupBY = pooledDf[self.indepVars].groupby('route')
            meanBY = groupBY.mean()
            XDemeaned = pooledDf[self.indepVars] - meanBY
            for id in self.idVar:
                tempData = XDemeaned[pooledDf[pooledDf.columns[0]] == id]







            self.bse = None
            self.cov_params = None
            self.params = None
            self.feParams = None
            self.pvalues = None
            self.betas = None
            self.Sigma = None
            self.SigmaRobust = None



            robustFE = FE.get_robustcov_results(cov_type='HC0', use_t=True)
            return robustFE
        return FE

    def betasREBalanced(self):
        #first estimate the pooled model and obtain the residuals
        pooledDF = self.panel.to_frame().reset_index()
        pooledModel = sm.ols(formula = self.formula, data = pooledDF).fit()
        sigmaSqV = pooledModel.ssr/pooledModel.df_resid
        sigmaSqC = 0
        tTime = self.timeVar[:-1]
        sTime = self.timeVar[1:]
        for id in self.idVar:
            for t in tTime:
                for s in sTime:
                    if t == s:
                        continue
                    sigmaSqC += pooledModel.resid[id].loc[t]*pooledModel.resid[id].loc[s]
        sigmaSqC /=(((pooledDF.shape[0]*(len(self.timeVar) - 1))/2) - pooledModel.df_model)
        sigmaSqU = sigmaSqV - sigmaSqC
        tempSigmaSqC = np.matrix(np.ones((len(self.timeVar), len(self.timeVar)))*sigmaSqC)
        tempSigmaSqU = np.matrix(np.identity(len(self.timeVar))*sigmaSqU)
        self.omega = tempSigmaSqC + tempSigmaSqU
        tempA = np.matrix(np.zeros((len(self.indepVars), len(self.indepVars))))
        tempXOY = np.matrix(np.zeros((len(self.indepVars), 1)))
        for id in self.idVar:
            tempData = np.matrix(self.panel.major_xs(id)[self.indepVars].as_matrix())
            tempY = np.matrix(self.panel.major_xs(id)[self.depVar].as_matrix())
            tempA += tempData.T * omega.inv * tempData
            tempXOY += tempData.T * omega.inv * tempY
        self.betas = tempA.inv * tempXOY
        self.Sigma = tempA.inv/len(self.idVar)
        if self.robust:
            tempB = np.matrix(np.zeros((len(self.indepVars), len(self.indepVars))))
            for id in self.idVar:
                tempData = np.matrix(self.panel.major_xs(id)[self.indepVars].as_matrix())
                tempY = np.matrix(self.panel.major_xs(id)[self.depVar].as_matrix())
                tempResid = tempY - tempData*self.betas
                tempB += tempData.T * self.omega.inv * tempResid * tempResid.T * self.omega.inv * tempData
            self.SigmaRobust = (tempA.inv * tempB * tempA.inv) / len(self.idVar)






class Results:
    def __init__(self):



        ###attributes to be defined later:
        self.bse = None
        self.cov_params = None
        self.params = None
        self.feParams = None
        self.pvalues = None

    def conf_int(self):
        pass

    def predict(self):
        pass

    def summary(self):
        pass