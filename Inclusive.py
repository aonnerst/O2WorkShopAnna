import numpy as np
import ROOT
import pickle #import 2.76 data
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt

import scipy
from scipy import interpolate

import sys
sys.path.append("JPyPlotRatio");
#sys.path.append("/home/jasper/Asiakirjat/projects/JPyPlotRatio");


import JPyPlotRatio

f = ROOT.TFile("AnalysisResults.root","read")

dataTypePlotParams = [
        {'plotType':'data','color':'C2','fmt':'s','xshift':0.0,'markersize':3.0},
        
]

# define panel/xaxis limits/titles
ny = 1;
nx = 2;
xlimits = [(0.,60.)];
ylimits = [(-2.7*1e-7,3.4*1e-7),(-1,800)];
xtitle = ["$p_{T}$","$\\eta$"];
ytitle = ["$1/N_{evt} dN_{ch}/dp_T$","$dN_{ch}/d\\eta$"];
Nevt = 7640.0;


plot = JPyPlotRatio.JPyPlotRatio(panels=(ny,nx),panelsize=(4,5),disableRatio=[0],
        panelPrivateScale=[1],
        panelPrivateRowBounds={1:ylimits[1]},
        legendPanel=0,legendLoc=(0.18,0.72),legendSize=7,ylabel=ytitle[0],ylabelRight=ytitle[1], xlabel={0:xtitle[0],1:xtitle[1]});


##plot.EnableLatex(True);

plot.GetAxes(1).yaxis.tick_right()

hist_pt = f.Get("twoparcorcombexample/hChargedPt")
hist_eta = f.Get("twoparcorcombexample/hDeltaEta")

hist_pt.Scale(1./Nevt,"width")
hist_eta.Scale(1.,"width")

p = plot.Add(0,hist_pt,**dataTypePlotParams[0]);
p = plot.Add(1,hist_eta,**dataTypePlotParams[0]);

plot.GetAxes(0).set_yscale("log")
for k,a in enumerate(plot.ax[:,0]):
        a.set_yscale("log");
        a.set_yticks([(10**i) for i in range(0,5)]);

plot.Plot();
plot.Save("figs/Inclusive.pdf");