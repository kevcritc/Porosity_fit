#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:05:16 2022

@author: phykc
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit


class Sio2Aunr:
    def __init__(self, filename, rows=7):
        self.filename=filename
        # Load data from excelfile
        self.df=pd.read_excel(self.filename)
        self.rows=rows
        self.water=np.array([100,80,70,60, 50, 40, 30])
        self.em=np.array([1.768,1.844, 1.885, 1.923, 1.965, 2.005, 2.047])
        self.glycerol=100-self.water
        self.ngly=1.46716
        self.nwater=1.333
        self.esio2=2.1
        self.eh20=1.77
        self.egly=1.9
        self.refractive_index=(self.water*self.nwater+self.glycerol*self.ngly)/100
    def wavelengthtoindex(self,WL):
        xdata=self.df.iloc[:,0].to_list()
        m=(xdata[-1]-xdata[0])/len(xdata)
        x=int((WL-xdata[0])/m)
        print(x)
        return x
    def guassfunc(self,x,xc,w,A):
        return A*np.exp(-0.5*(x-xc)**2/w**2)
    # convert Re(e1) into wavelength
    def Landa(self,e1):
        return (e1-43.718)/-0.0847
    # convert wavelength into Re(e1) 
    def e1(self,landa):
        return -0.0847*landa+43.718
    # determine L3 from aspect ratio
    def L3(self,ar):
        return (1+ar)**-1.6
    # determine 'x' from L3
    def x(self,ar):
        return (1-self.L3(ar))/self.L3(ar)
    # plot the uv-vis data in spreadsheet
    def show_uvvis(self):
        for n in range(1,self.rows+1):
            plt.plot(self.df.iloc[:,0],self.df.iloc[:,n])
        plt.show()
    def f1(self,emed):
        return -1466.7*emed**5 + 14289*emed**4 - 55668*emed**3 + 216819*emed**2 - 211067*emed + 82168
 
    def loadpeaks(self, df, dataset):
        self.df_peaks=df
        ind_val=dataset-1
        colindex=3*ind_val+3
        self.peakfitlist=self.df_peaks.iloc[:,colindex].tolist()
        self.error=self.df_peaks.iloc[:,colindex+1].tolist()
        
    def findpeaks(self, aboveWL=650):
        self.frompoint=self.wavelengthtoindex(aboveWL)
        # Find rough peak position using peakfind
        self.peaks_list=[]
        for n in range(1,self.rows+1):
            peaks=find_peaks(self.df.iloc[:,n])
            maxima=self.df.iloc[peaks[0],0].to_list()
            self.peaks_list.append(maxima)
        # Using the rough peak positions fit a guassian to the LSPR to get an accuate peak pos.
        self.peakfitlist=[]
        self.error=[]
        for n in range(1,self.rows+1):
            xdata=self.df.iloc[self.frompoint:,0].to_list()
            ydata=self.df.iloc[self.frompoint:,n].to_list()
        
            popt, pcov =curve_fit(self.guassfunc,xdata,ydata, p0=[self.peaks_list[n-1][-1],100,0.2])
            for i in range(0,len(popt),3):
                try:
                  self.error.append(np.absolute(pcov[i][i])**0.5)
                except:
                  self.error.append( 0.00 )
            self.peakfitlist.append(popt[0])
    def plotpeaks(self):
        # Plot data versus n
        plt.scatter(self.refractive_index,self.peakfitlist)
        plt.xlabel('n, refractive index')
        plt.ylabel('LSPR Peak maxima /nm')
        plt.show()
    def mod_LSPR(self,em,ar,f):
        x_val=self.x(ar)
        emod=self.esio2*(1-f)+em*f
        E1=-1*emod*x_val
        WLs=self.Landa(E1)
        return WLs
    def fit_model(self):
        print('list of peak positions: ',self.peakfitlist)
        # Next fit data to modelled effective dielectic constant data with AR and f as fit aparameters
        self.LSPR_peaks=np.array(self.peakfitlist)
        # Input the medium dielect constants
        
        self.popt2, self.pcov2=curve_fit(self.mod_LSPR,self.em,self.LSPR_peaks, p0=[4,0.8])
        self.error1=[]
        # Determine the error of the fits
        for i in range(0,len(self.popt2)):
            try:
              self.error1.append(np.absolute(self.pcov2[i][i])**0.5)
            except:
              self.error1.append( 0.00 )
        print('Best fit, AR= ',self.popt2[0],' and f= ',self.popt2[1])
        print('errors are:')
        print(self.error1)
        # Plot the data
    def plotfit(self, COl, labelname):
        plt.errorbar(self.em,self.LSPR_peaks,yerr=self.error, marker='o',linestyle='', capsize=3, label=labelname)
        plt.plot(self.em,self.mod_LSPR(self.em,*self.popt2),color=COl)
        plt.xlabel(r'e$_{m}$')
        plt.ylabel('LSPR maxima /nm')
        

if __name__=="__main__":
    filename1='porosity -no CTAB.xlsx'
    filename2='porosity -5 mM CTAB.xlsx'
    filename3='porosity -10 mM CTAB.xlsx'
    filename4='porosity -15 mM CTAB.xlsx'
    
    filepeaks='All LSPR peak positions 21-7-22.xlsx'
    df=pd.read_excel(filepeaks)
    
    porSiO21=Sio2Aunr(filename1, rows=7)
    porSiO21.show_uvvis()
    porSiO21.findpeaks()
    porSiO21.loadpeaks(df, dataset=1)
    porSiO21.fit_model()
    
    porSiO2=Sio2Aunr(filename2, rows=7)
    porSiO2.show_uvvis()
    porSiO2.findpeaks()
    porSiO2.loadpeaks(df, dataset=2)
    porSiO2.fit_model()
    
    por3=Sio2Aunr(filename3, rows=7)
    por3.show_uvvis()
    por3.findpeaks()
    por3.loadpeaks(df, dataset=3)
    por3.fit_model()
    
    por4=Sio2Aunr(filename4, rows=7)
    por4.show_uvvis()
    por4.findpeaks()
    por4.loadpeaks(df, dataset=4)
    por4.fit_model()
    
    porSiO21.plotfit('blue','0 mM CTAB')
    porSiO2.plotfit('red', '5 mM CTAB')
    por3.plotfit('green','10 mM CTAB')
    por4.plotfit('orange','15 mM CTAB')
    plt.legend(loc=4)
    plt.show()



