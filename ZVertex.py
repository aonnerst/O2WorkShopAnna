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
        {'plotType':'data','color':'C1','fmt':'s','markersize':2.0},
]


# define panel/xaxis limits/titles
ny = 1;
nx = 1;
xlimits = [(0.,60.)];
xtitle = ["Z-vertex (cm)"];
ytitle = ["$dN_{evt}/dz$"];


plot = JPyPlotRatio.JPyPlotRatio(panels=(ny,nx),panelSize=(4,5),disableRatio=[0],
	legendPanel=0,legendLoc=(0.18,0.72),legendSize=7,ylabel=ytitle,xlabel=xtitle);

##plot.EnableLatex(True);

hist_zvtx =f.Get("twoparcorcombexample/hZvertex")
hist_zvtx.Print()

nBins = 200
Nevt = hist_zvtx.GetEntries()
print(Nevt)
hist_zvtx.Scale(1.,"width")
p = plot.Add(0,hist_zvtx,**dataTypePlotParams[0]);

plot.Plot();
plot.Save("figs/ZVertex.pdf");
#plot.Show();	

f.Close();