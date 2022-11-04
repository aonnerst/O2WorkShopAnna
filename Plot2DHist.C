//#include "include/const_2pc.h"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>

void Plot2DHist(){
	TFile* fin = new TFile("AnalysisResults.root", "read");
	TDirectoryFile* tdf = (TDirectoryFile*)fin->Get("twoparcorcombexample");
	TH2F* hDeltaPhiDeltaEta = (TH2F*) tdf->Get("hDeltaPhiDeltaEta");
	TH1F* hTrigg = (TH1F*)tdf->Get("hTrigg");

	TCanvas* c = new TCanvas("c","c",600,600);

	int Ntriggers=hTrigg->GetEntries();

	hDeltaPhiDeltaEta->Scale(1.0/Ntriggers,"width");

	hDeltaPhiDeltaEta->Draw("surf1");

	c-> Update();
	c->Show();
}