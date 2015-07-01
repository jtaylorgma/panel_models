########################################################################################################################
#Program: testPanel.py
#Project: Panel Models
#Author: Josh Taylor
#Last Edited: 6/30/15
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
from numpy.linalg import inv
import scipy as sp
from math import sqrt
import re


#path to this file: C:\Users\jtaylor\Projects\Prais-Winsten\Python\Code\panel_models\testPanel.py


class PanelModel:



    def __init__(self, formula, effects = "random", intercept = True, time_fe = False, entity_fe = False, robust = False, data = None):
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
        self.pooledDF = self.panel.to_frame(filter_observations = False).reset_index()
        self.balanced = True
        self.timeVals = self.panel.minor_axis
        self.idVals = self.panel.major_axis
        self.indepVars = re.findall(r"\w+",formula)[1:]
        self.depVar = re.findall(r"\w+",formula)[0]
        self.intercept = intercept
        if self.intercept:
            _needIntercept = True
            for v in self.indepVars:
                if len(pd.unique(self.panel[v].values.ravel())) == 1:
                    _needIntercept = False
            if _needIntercept:
                self.panel['constant'] = 1
                self.indepVars.append('constant')
        self.idVar = self.panel.to_frame().reset_index().columns[0]
        self.timeVar = self.panel.to_frame().reset_index().columns[1]

        if time_fe:
            if self.effects == "fixed":
                self.formula = self.formula + " + C(%s)" %self.timeVar
            else:
                for s in self.timeVals:
                    self.pooledDF['t%s' %s] = self.pooledDF[self.timeVar] == s
                    if s == self.timeVals[0]:
                        continue
                    self.formula += " + t%s" %s
                    self.indepVars.append('t%s' %s)
                _tempPanel= self.pooledDF.set_index([self.idVar, self.timeVar])
                self.panel = _tempPanel.to_panel()

        if entity_fe and effects == "random":
            self.formula = self.formula + " + C(%s)" %self.idVar

        ### attributes to be defined later:
        self.params = None
        self.bse = None
        self.params = None
        self.feParams = None
        self.pvalues = None
        self.cov_params = None
        self.cov_params_robust = None
        self.omega = None
        self.tvalues = None
        self.conf_int = None
        self.resid = None
        self.fittedvalues = None
        self.rsquared = None
        self.ssr = None
        self.sst = None
        #TODO: RE unbalanced

        #TODO: Rethink how you want the fixed effects defined

    def balanceChecker(self):
        message = "This panel is balanced"
        for id in self.idVals:
            for t in self.timeVals:
                if t not in self.panel.major_xs(id).index: #checking if each id has all times
                    message = "This panel is unbalanced"
                    self.balanced = False
                    break
        print message

    def fit(self):
        if self.effects == "fixed":
            self.fixedEffects()
            return self
        elif self.effects == 'random':
            self.randomEffects()
            return self

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
        formulaFE = self.formula + " + C(%s)" %self.idVals
        FE = sm.ols(formula=formulaFE, data = pooledDf).fit()

        self.params = FE.params.get(self.indepVars)
        self.bse = FE.bse.get(self.indepVars)
        self.feParams = FE.params.drop(self.indepVars)
        self.feParamsPValues = FE.pvalues.drop(self.indepVars)
        self.pvalues = FE.pvalues.get(self.indepVars)
        self.cov_params = FE.cov_params
        self.tvalues = FE.tvalues
        self.conf_int = FE.conf_int()
        self.resid = FE.resid
        self.fittedvalues = FE.fittedvalues
        self.mse = FE.mse_resid
        self.ssr = FE.ssr
        self.rsquared = FE.rsquared



        if self.robust:
            pooledDf['resid'] = FE.resid
            groupBY = pooledDf[self.indepVars].groupby('route')
            meanBY = groupBY.mean()

            for id in self.idVals:
                XDemeaned = pd.DataFrame()
                tempData = pooledDf[pooledDf[self.idVar] == id].iloc[:,:-1]
                sumXDemeaned = pd.DataFrame(index=len(self.timeVals), columns=len(pooledDf.columns) - 2)
                tempSumXDemeaned = pd.DataFrame(index=len(self.timeVals), columns=len(pooledDf.columns) - 2)
                sumXDemeaned = sumXDemeaned.fillna(0)
                tempSumXDemeaned = tempSumXDemeaned.fillna(0)
                for i in range(tempData.shape[0]):
                    dmTempData = tempData.iloc[i,1:] - meanBY.loc[id,:-1]
                    dmTempData[self.idVar] = id
                    XDemeaned = XDemeaned.append(dmTempData)
                tempSumXDemeaned += np.matrix(XDemeaned[XDemeaned[self.idVar] == id].as_matrix())
                tempResid = np.matrix(pooledDf[pooledDf[self.idVar] == id].iloc[:,-1].as_matrix())
                sumXDemeaned += tempSumXDemeaned.T * tempResid * tempResid.T * tempSumXDemeaned

            X = np.matrix(XDemeaned.as_matrix())
            self.cov_params_robust = inv(X.T * X) * sumXDemeaned * inv(X.T * X)
            self.bse = pd.Series()
            for i in range(len(self.indepVars)):
                self.bse[self.indepVars[i]] = sqrt(self.cov_params_robust[i,i])
            self.tvalues = self.params / self.bse
            _df = FE.nobs - FE.df_model - 1
            _dfVec = np.ones(len(self.indepVars)) * _df
            self.pvalues = sp.stats.tvalues.sf(self.tvalues, _dfVec)*2
            _ciLowerBound = self.params - 1.96*self.bse
            _ciUpperBound = self.params + 1.96*self.bse
            self.conf_int = pd.merge(_ciLowerBound.to_frame(), _ciUpperBound.to_frame(), left_index=True, right_index=True)





    def betasREBalanced(self):
        #first estimate the pooled model and obtain the residuals
        pooledModel = sm.ols(formula = self.formula, data = self.pooledDF).fit()
        sigmaSqV = pooledModel.ssr/pooledModel.df_resid
        sigmaSqC = 0
        tTime = self.timeVals.tolist()[:-1]
        idList = self.idVals.tolist()
        timeList = self.timeVals.tolist()
        for id in self.idVals:
            for t in tTime:
                tIndex = self.timeVals.tolist().index(t)
                sTime = self.timeVals.tolist()[(tIndex + 1):]
                for s in sTime:
                    sigmaSqC += pooledModel.resid[idList.index(id)*len(timeList) + timeList.index(t)]*pooledModel.resid[idList.index(id)*len(timeList) + timeList.index(s)]
        sigmaSqC /=(((self.pooledDF.shape[0]*(len(self.timeVals) - 1))/2) - pooledModel.df_model)
        sigmaSqU = sigmaSqV - sigmaSqC
        tempSigmaSqC = np.matrix(np.ones((len(self.timeVals), len(self.timeVals)))*sigmaSqC)
        tempSigmaSqU = np.matrix(np.identity(len(self.timeVals))*sigmaSqU)
        self.omega = tempSigmaSqC + tempSigmaSqU
        tempA = np.matrix(np.zeros((len(self.indepVars), len(self.indepVars))))
        tempXOY = np.matrix(np.zeros((len(self.indepVars), 1)))
        for id in self.idVals:
            tempData = np.matrix(self.panel.major_xs(id)[self.indepVars].as_matrix())
            tempY = np.matrix(self.panel.major_xs(id)[self.depVar].as_matrix())
            tempA += tempData.T * inv(self.omega) * tempData
            tempXOY += tempData.T * inv(self.omega) * tempY.T
        _params = (inv(tempA) * tempXOY).A1
        self.params = pd.Series(_params, self.indepVars)
        _cov_params = inv(tempA)
        self.cov_params = pd.DataFrame(_cov_params, self.indepVars, self.indepVars)
        self.bse = pd.Series()
        for i in range(len(self.indepVars)):
            self.bse[self.indepVars[i]] = sqrt(_cov_params.item((i,i)))
        self.tvalues = self.params/self.bse
        _df = np.ones(len(self.indepVars))*(self.pooledDF.shape[0] - len(self.indepVars) - 2)
        self.pvalues = sp.stats.t.sf(self.tvalues, _df)*2
        self.fittedvalues = np.dot(self.pooledDF[self.indepVars].as_matrix(), (np.array(self.params)).T)
        self.resid = self.pooledDF[self.depVar] - self.fittedvalues
        self.mse = np.sum(self.resid**2)/(self.pooledDF.shape[0] - len(self.indepVars) - 2)
        self.rmse = sqrt(self.mse)
        self.ssr = np.sum(self.resid**2)
        _yBar = np.mean(self.pooledDF[self.depVar])
        self.sst = np.sum(self.pooledDF[self.depVar] - _yBar)
        self.rsquared = 1 - self.ssr/self.sst

        if self.robust:
            tempB = np.matrix(np.zeros((len(self.indepVars), len(self.indepVars))))
            for id in self.idVals:
                tempData = np.matrix(self.panel.major_xs(id)[self.indepVars].as_matrix())
                tempY = np.matrix(self.panel.major_xs(id)[self.depVar].as_matrix())
                tempResid = tempY - tempData*self.params
                tempB += tempData.T * inv(self.omega) * tempResid * tempResid.T * inv(self.omega) * tempData
            _cov_params_robust = (inv(tempA) * tempB * inv(tempA)) / len(self.idVals)
            self.cov_params_robust = pd.DataFrame(_cov_params_robust, self.indepVars, self.indepVars)
            for i in range(len(self.indepVars)):
                self.bse[self.indepVars[i]] = sqrt(self.cov_params_robust[i,i])
            self.tvalues = self.params/self.bse
            self.pvalues = sp.stats.tvalues.sf(self.tvalues, _df)*2

        _ciLowerBound = self.params - 1.96*self.bse
        _ciUpperBound = self.params + 1.96*self.bse
        self.conf_int = pd.merge(_ciLowerBound.to_frame(), _ciUpperBound.to_frame(), left_index=True, right_index=True)


    def summary(self):
        pass

        #This is going to have to wait


        """
        trunParams = self.params
        trunBSE = self.bse
        trunTValue = self.tvalues
        for p in range(len(trunParams)):
            trunParams[p] = float(str(trunParams[p])[:6]) if len(str(trunParams[p])) > 6 else trunParams[p]
            trunBSE[p] = float(str(trunBSE[p])[:6]) if len(str(trunBSE[p])) > 6 else trunBSE[p]
            trunTValue = float(str(trunTValue[p])[:4]) if len(str(trunTValue[p])) > 4 else trunBSE[p]



        print
        """



cornwell = pd.read_csv("C:\Users\jtaylor\Downloads\cornwell.csv")
cornwell['constant'] = 1
cornwell = cornwell.set_index(['county', 'year'])
cornwellPanel = cornwell.to_panel()
cornwellRE1 = PanelModel(formula='crmrte ~ polpc + urban + prbpris', data= cornwellPanel).fit()
cornwellRE2 = PanelModel(formula='crmrte ~ polpc + urban + prbpris', data= cornwellPanel, time_fe=True).fit()























