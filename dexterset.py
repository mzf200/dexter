#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:47:00 2019

@author: maxf
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class dexterset:
      
      def __init__(self,filepath,commonvar = 0,runs = 'detect'):
            
            self.__unique_array = None
            self.__col_names = ['File Number','Variable', 'N [OD]', 'Nx','Nz',\
                                'Npx','Tx','Tz','PSD', 'OD','n_pk','sigma_x',\
                                'sigma_z,','x','z']
            
            self.__df = pd.read_csv(filepath,header=None,\
                                    names=self.__col_names)
            self.__filenums = self.__df[self.__col_names[0]].values.tolist()
            self.__vars = self.__df[self.__col_names[1]].values.tolist()
            self.__nods = self.__df[self.__col_names[2]].values.tolist()
            self.__nxs = self.__df[self.__col_names[3]].values.tolist()
            self.__nzs = self.__df[self.__col_names[4]].values.tolist()
            self.__npxs = self.__df[self.__col_names[5]].values.tolist()
            self.__Txs = self.__df[self.__col_names[6]].values.tolist()
            self.__Tzs = self.__df[self.__col_names[7]].values.tolist()
            self.__psds = self.__df[self.__col_names[8]].values.tolist()
            self.__ods = self.__df[self.__col_names[9]].values.tolist()
            self.__npks = self.__df[self.__col_names[10]].values.tolist()
            self.__sigxs = self.__df[self.__col_names[11]].values.tolist()
            self.__sigzs = self.__df[self.__col_names[12]].values.tolist()
            self.__xs = self.__df[self.__col_names[13]].values.tolist()
            self.__zs = self.__df[self.__col_names[14]].values.tolist()
            self.__data_array = np.zeros((len(self.__vars),len(self.__col_names)))
            self.__data_array[:,0] = self.__filenums
            self.__data_array[:,1] = self.__vars
            self.__data_array[:,2] = self.__nods
            self.__data_array[:,3] = self.__nxs
            self.__data_array[:,4] = self.__nzs
            self.__data_array[:,5] = self.__npxs
            self.__data_array[:,6] = self.__Txs
            self.__data_array[:,7] = self.__Tzs
            self.__data_array[:,8] = self.__psds
            self.__data_array[:,9] = self.__ods
            self.__data_array[:,10] = self.__npks
            self.__data_array[:,11] = self.__sigxs
            self.__data_array[:,12] = self.__sigzs
            self.__data_array[:,13] = self.__xs
            self.__data_array[:,14] = self.__zs
            
            
# =============================================================================
#             if runs == 'detect':
#                   
#                   if type(commonvar) == int:
#                         intcount = self.__vars.count(commonvar)
#                         floatcount = self.__vars.count(float(commonvar))
#                         runs = max([intcount,floatcount])
#                   
#                   else:
#                         runs = self.__vars.count(commonvar)
#                  
#             self.__runs = runs
#             
#             self.__runlist = [] #[[] for j in range(self.__runs)]
#             
# =============================================================================
            
            runslice = []
            runlist = []
            for j in range(len(self.__vars)):
                  if self.__vars[j] in runlist:
                        self.__runlist.append(runslice)
                        runslice = [j]
                        runlist = [self.__vars[j]]
                        
                  else:
                        runslice.append(j)
                        runlist.append(self.__vars[j])
            
            self.__runlist = runlist
            self.__runs = len(runlist)
         
#==============================================================================                                
#==============================================================================
                        
      def data_sort(self):
            
            sort_vars = np.sort(self.__vars)
            unique_vars = np.unique(sort_vars)
            self.__unique_array = np.zeros((unique_vars.shape[0],self.__data_array.shape[1]))
            self.__std_array = np.zeros((unique_vars.shape[0],self.__data_array.shape[1]))
            xval = unique_vars[0]
            temp = [[] for j in range(self.__std_array.shape[1]-1)]
            
            for j in range(len(self.__vars)):
                  if self.__vars[j] == xval:
                        for k in range(1,self.__data_array.shape[1]):
                              temp[k-1].append(self.__data_array[j,k])
                              
            mean_list = [np.mean(k) for k in temp]
            std_list = [np.std(k)/len(k) for k in temp]
            
            self.__unique_array[0,0] = xval
            self.__unique_array[0,1:] = mean_list
            self.__std_array[0,0] = xval
            self.__std_array[0,1:] = std_list            
            
            
            for j in range(1,unique_vars.shape[0]):
                  xval = unique_vars[j]
                  temp = [[] for j in range(self.__std_array.shape[1]-1)]
                  for k in range(len(self.__vars)):
                        if self.__vars[k] == xval:
                              for m in range(1,self.__data_array.shape[1]):
                                    temp[m-1].append(self.__data_array[k,m])
                                    
                  mean_list = [np.mean(k) for k in temp]
                  std_list = [np.std(k)/len(k) for k in temp]
                  self.__unique_array[j,0] = xval
                  self.__unique_array[j,1:] = mean_list
                  self.__std_array[j,0] = xval
                  self.__std_array[j,1:] = std_list
                                          
            return self.__unique_array, self.__std_array
      
#==============================================================================
#==============================================================================
      
      def data_array(self):
            return self.__data_array
      
#==============================================================================
#==============================================================================
      
      def dexterplot(self,xvariable,yvariable,std_err=False,pointstyle = 'bx'):
            if self.__unique_array ==  None:
                  self.data_sort()
            
            for j,k in enumerate(self.__col_names):
                  if k == xvariable:
                        xindex = j
                  if k == yvariable:
                        yindex = j
            
            plt.figure()
            if std_err == True:
                  plt.errorbar(self.__unique_array[:,xindex],
                               self.__unique_array[:,yindex],fmt=pointstyle,
                               yerr = self.__std_array[yindex,:])
                  
            else:
                  plt.plot(self.__unique_array[:,xindex],
                           self.__unique_array[:,yindex],pointstyle)
                  
#==============================================================================
#==============================================================================
                  
      def dexter_multiplot(self,xvariable,yvariable,std_err=False,runs = 'all',
                           pointstyle = 'x'):
            
            for j,k in enumerate(self.__col_names):
                  if k == xvariable:
                        xindex = j
                  if k == yvariable:
                        yindex = j
            
            plt.figure()
            
            for j in range(self.__runs):
                  xvar = [self.__data_array[k,xindex] for k in self.__runlist[j]]
                  yvar = [self.__data_array[k,yindex] for k in self.__runlist[j]]
            
                  if std_err == True:
                        error = [self.__std_array[k,yindex-1] for k in self.__runlist[j]]
                        
                        plt.errorbar(xvar,yvar,yerr = error,fmt = pointstyle)
                        
                  else:
                        
                        plt.plot(xvar,yvar,pointstyle)
                  
            plt.show()
            
                  
      
