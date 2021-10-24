//  Created by Miles Cochran-Branson on 7/19/21.
//
//  Calculation of final production of quarkonium pairs
//  Uses values from the pythia_simulations folder calculated from PYTHIA
//
//  This file analyzes data from the 2018 PbPb run at CMS
//  This run has an integrated luminosity of 1.7nb^-1 and sqrt(s) = 5.02 TeV
//  Specifically, we look at the production of Upsilon and J/psi particles
//  as well as pairs of particles being produced. 
//  This script makes a histogram showing the values of predicted production
//  of these particles and pairs of particles along with uncertainties. 
//
//  REQUIRES:
//  A glauber model and .root file with Ncoll data, data_values.h header. 
//  Found in data_values.h: 
//  Upsilon and J/psi probability of production,
//  PbPb cross-section from Glauber model,
//  efficiency and acceptance measurements. 

//Calculated Values here:
#include <iostream>
#include <TChain.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include "data_values.hpp"

using namespace std; 

void jpsiupdouble() {

    //Instance of class
    DataVals vals;
    
    //Get values from .root file
    TChain in_chain("glauberfinal70mb10M");
    in_chain.Add("glauberfinal70mb10M*.root");
    
    //Set number of Bins
    int numBins = 2400;
    
    //Define hist. to hold Ncoll
    TH1F* numColl = new TH1F("numCollisions","Number Upsilons Produced per N_{coll}; N_{coll}; counts", numBins, -0.5, 2399.5);
    
    float b, NColl, Npart;
    in_chain.SetBranchAddress("NumCollisions", &NColl);
    
    //Fill histo with data!
    for (int irow = 0; irow<in_chain.GetEntries(); ++irow) {
        in_chain.GetEntry(irow);
        
        numColl->Fill(NColl);
        
    }

    //Definition of Probabilities 
    double JpsiProb = vals.JpsiProb();
    double UpsilonProb = vals.UpsProb();
    
    double prob[3];
    prob[1] = JpsiProb; //prob of j/psi production
    prob[2] = UpsilonProb; //prob of upsilon production
    prob[0] = 1 - (prob[1] + prob[2]); //prob of no special events

    //define acceptance and efficieny 
    bool useAcceptEff = vals.useAcceptEff;
    double JpsiA_E = vals.JpsiAE();
    double UpsilonA_E = vals.UpsAE();

    double accepteff[6];
    if (useAcceptEff == true) {
        accepteff[0] = 1, accepteff[1] = JpsiA_E, accepteff[2] = UpsilonA_E, accepteff[3] = JpsiA_E * JpsiA_E;
        accepteff[4] = JpsiA_E * UpsilonA_E; accepteff[5] = UpsilonA_E * UpsilonA_E;
    }
    else {
        for (int i(0); i < 6; ++i) accepteff[i] = 1;
    }

    const char *title;
    if (useAcceptEff == true) title = "#scale[0.75]{#splitline{Production of J/#psi and #varUpsilon}{with Acceptance and Efficiency}}; ;Events";
    else title = "#scale[0.75]{#splitline{Production of J/#psi and #varUpsilon}{without Acceptance and Efficiency}}; ;Events";

    TH1F* ParticleProductionHist = new TH1F("ParticleProductionHist", title, 6, -0.5, 5.5);
    
    double production[6] = {0.0};
    double scale;
    double NEvennts = 0;
    
    for (int i(0); i < numBins; ++i) {
        
        scale = numColl->GetBinContent(i);

        NEvennts += scale * i;
        
        production[0] += pow(prob[0], i) * scale; //no special particle production
        production[1] += i * pow(prob[0], i - 1) * prob[1] * scale; //single j/psi
        production[2] += i * pow(prob[0], i - 1) * prob[2] * scale; //single upsilon
        production[3] += (i / 2) * (i - 1) * pow(prob[0], (i - 2)) * pow(prob[1], 2) * scale; //double j/psi
        production[4] += i * (i - 1) * pow(prob[0], (i - 2)) * prob[1] * prob[2] * scale; //j/psi + upsilon
        production[5] += (i / 2) * (i - 1) * pow(prob[0], (i - 2)) * pow(prob[2], 2) * scale; //double upsilon
    }
    
    for (int j(0); j < 6; ++j) {
        production[j] = production[j] * accepteff[j]; //add acceptance and efficiency into calculation
        ParticleProductionHist->SetBinContent(j + 1, production[j]);
    }
    
    //re-scale to approprieate number of events
    ParticleProductionHist->Scale(1.0/ParticleProductionHist->Integral());
    ParticleProductionHist->Scale(vals.NumEvents()); //scale from calculated PbPbsigma and L_int
    
    
    //Draw Commands, etc
    gStyle->SetOptStat(0);
    gStyle->SetLegendTextSize(0.04);
    
    TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
    c1->SetLogy();
    
    //labels for histogram 
    const char *label[] = {"none", "J/#psi", "#varUpsilon", "J/#psi + J/#psi", "J/#psi + #varUpsilon", "#varUpsilon + #varUpsilon"};
    
    for (int k(0); k < 6; ++k) {
        ParticleProductionHist->GetXaxis()->SetBinLabel(k + 1, label[k]);
    }
    
    //garph settings
    ParticleProductionHist->GetXaxis()->SetLabelSize(0.05);
    ParticleProductionHist->GetYaxis()->SetTitleOffset(1.4); //adjust placement of y title
    ParticleProductionHist->SetLineColor(kViolet-6); 
    ParticleProductionHist->SetLineWidth(5);
    ParticleProductionHist->SetFillColor(kMagenta-8); //fill with slightly darker color
    ParticleProductionHist->Draw("hist");
    
    //legend settings
    TLegend* leg = new TLegend(0.75, 2e+5, 5.5, 2.5e+10, "", "");
    leg->AddEntry(ParticleProductionHist, "Multinomial Model");
    leg->AddEntry((TObject*)0, "Total Events = 1.28 #times 10^{10}", "");
    leg->AddEntry((TObject*)0, Form("%s #approx %.0f #pm %.0f", label[0], ParticleProductionHist->GetBinContent(1), ParticleProductionHist->GetBinError(1)), "");
    for (int i(1); i < 6; ++i) {
        leg->AddEntry((TObject*)0, Form("%s #approx %.3f #pm %.3f", label[i], ParticleProductionHist->GetBinContent(i+1), ParticleProductionHist->GetBinError(i+1)), "");
    }
    leg->Draw();   

    int whichProb = vals.whichProb;

    //latex text
    TLatex TeX;
    const char* probabilityTitle;
    if (whichProb == 1) probabilityTitle = "#scale[0.8]{#splitline{Pr(#varUpsilon) = 3.10 #times 10^{-7}}{Pr(J/#psi) = 3.93 #times 10^{-5}}}";
    if (whichProb == 2) probabilityTitle = "#scale[0.8]{#splitline{Pr(#varUpsilon) = 1.34 #times 10^{-7}}{Pr(J/#psi) = 3.91 #times 10^{-5}}}";
    if (whichProb == 3) probabilityTitle = "#scale[0.8]{#splitline{Pr(#varUpsilon) = 2.26 #times 10^{-5}}{Pr(J/#psi) = 3.28 #times 10^{-4}}}";
    else return;
    TeX.DrawLatex(2.5, 1e+7, probabilityTitle);
    TeX.DrawLatex(2.5, 5e+8, "#scale[0.8]{#int L = 1.7nb^{-1}}");
    TeX.DrawLatex(2.5, 5e+9, "#scale[0.8]{#sqrt{s_{NN}} = 5.02 TeV}");
    TeX.DrawLatex(2.5, 5e+8, "#scale[1.5]{#splitline{#varUpsilon + J/#psi Production}{with A #times E}}");

    //labels for histogram 
    const char *label1[] = {"none", "J/\\psi", "\\Upsilon", "J/\\psi + J/\\psi", "J/\\psi + \\Upsilon", "\\Upsilon + \\Upsilon"};
    cout << vals.NumEvents() << endl;
    cout << " " << endl;
    for (int i(0); i < 6; ++i) {
        cout << Form("$%s$ & $%.0f \\pm %.0f$ \\\\ ", label1[i], ParticleProductionHist->GetBinContent(i+1), ParticleProductionHist->GetBinError(i+1)) << endl;
    }

}

int main() {
    jpsiupdouble();
    return 0;
}