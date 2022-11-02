import numpy as np
import ROOT
import pickle #import 2.76 data
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
import math

import scipy
from scipy import interpolate

import sys
sys.path.append("JPyPlotRatio");
#sys.path.append("/home/jasper/Asiakirjat/projects/JPyPlotRatio");

from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F, TF1

import JPyPlotRatio

f = ROOT.TFile("AnalysisResults.root","read")

dataTypePlotParams = [
        {'plotType':'data','color':'C3','fmt':'s','xshift':0.0,'markersize':5.0},
        {'plotType':'data','color':'C1' ,'fmt':'s','xshift':-1.0,'markersize':5.0},
        {'plotType':'data','color':'C2','fmt':'o','xshift':-0.5,'fillstyle':'none','markersize':5.0},
        {'plotType':'data','color':'C3' ,'fmt':'o','xshift':0.0,'markersize':5.0},
        {'plotType':'data','color':'C4' ,'fmt':'D','xshift':0.5,'fillstyle':'none','markersize':5.0},
        {'plotType':'data','color':'C5','fmt':'D','xshift':1.0,'markersize':5.0},
        {'plotType':'data','color':'C6' ,'fmt':'8','xshift':1.5,'fillstyle':'none','markersize':5.0},
        {'plotType':'data','color':'C7','fmt':'8','xshift':2.0,'markersize':5.0},
        {'plotType':'data','color':'C8' ,'fmt':'h','xshift':2.5,'fillstyle':'none','markersize':5.0},
        {'plotType':'data','color':'C9' ,'fmt':'h','xshift':0.5,'markersize':5.0},
        {'plotType':'data','color':'blue','fmt':'s','fillstyle':'none','markersize':5.0},
        {'plotType':'data','color':'green' ,'fmt':'s','markersize':5.0},
        {'plotType':'data','color':'plum','fmt':'o','fillstyle':'none','markersize':5.0},
        
]

# define panel/xaxis limits/titles
ny = 1;
nx = 1;
xlimits = [(0.,60.)];
ylimits = [(-2.7*1e-7,3.4*1e-7),(-0.31,1.6)];
xtitle = ["Bin"];
ytitle = ["$<Cos(2*\\Delta\\phi)>$"];

plot = JPyPlotRatio.JPyPlotRatio(panels=(ny,nx),panelsize=(4,5),disableRatio=[0],

        legendPanel=0,legendLoc=(0.18,0.72),legendSize=7,ylabel=ytitle[0]);


##plot.EnableLatex(True);

hist_corr = f.Get("twoparcorcombexample/hTPcosDeltaPhi")
#Getting the avarage over event
average_corr = hist_corr.GetMean()
average_corr_err = hist_corr.GetMeanError()

vn = math.sqrt(abs(average_corr))
vn_err = 0.5*(1.0/(average_corr**0.5)*average_corr_err)
x = np.array(1.0)
y = np.array(vn)
x_err = np.array(0.0)
y_err = np.array(vn_err)
print(x,vn,vn_err)
text = "$v_2 = {:.03f}±{:.03f}$".format(vn,vn_err)

p = plot.Add(0,(x,y,y_err),**dataTypePlotParams[0]);

plot.GetPlot().text(0.17,0.75,text,fontsize=10);

plot.Plot();
plot.Save("figs/v2.pdf");