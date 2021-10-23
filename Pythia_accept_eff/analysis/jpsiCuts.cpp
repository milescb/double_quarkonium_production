// Author: Miles Cochran-Branson, 7/26/21
//
// cuts for j/psi data. Requires .root file created from PYTHIA simulation
// makes kinematic cuts to determine acceptance and efficiency 

void jpsiCuts() {

    //extract info from root files
    TChain chain("JpsiDat");
    chain.Add("JpsiData_final1*.root");
    
    int numBins = 100;

    //make histograms here 
    TH1F* JpsiHist[9];
    int k(0);
    JpsiHist[k++] = new TH1F("JpsiE", "J/#psi Energy; Energy [GeV]; Events", numBins, -0.5, 1999.5);
    JpsiHist[k++] = new TH1F("JpsipT", "J/#psi p_{T}; p_{T} [GeV/c]; Events", numBins, -0.5, 30);
    JpsiHist[k++] = new TH1F("JpsiRap", "J/#psi Rapidity; y; Events", numBins, -7.5, 7.5);
    JpsiHist[k++] = new TH1F("muonpE", "#mu^{+} Energy; Energy [GeV]; Events", numBins, -0.5, 1999.5);
    JpsiHist[k++] = new TH1F("muonppT", "#mu^{+} p_{T}; p_{T} [GeV/c]; Events", numBins, -0.5, 30);
    JpsiHist[k++] = new TH1F("muonpy", "#mu^{+} Rapidity; y; Events", numBins, -7.5, 7.5);
    JpsiHist[k++] = new TH1F("muonmE", "#mu^{-} Energy; Energy [GeV]; Events", numBins, -0.5, 1999.5);
    JpsiHist[k++] = new TH1F("muonmpT", "#mu^{-} p_{T}; p_{T} [GeV/c]; Events", numBins, -0.5, 30);
    JpsiHist[k++] = new TH1F("muonmy", "#mu^{-} Rapidity; y; Events", numBins, -7.5, 7.5);
    
    //Set Branch Addresses 
    float entry[9];
    chain.SetBranchAddress("JpsiE", &entry[0]);
    chain.SetBranchAddress("JpsipT", &entry[1]);
    chain.SetBranchAddress("Jpsiy", &entry[2]);
    chain.SetBranchAddress("mupE", &entry[3]);
    chain.SetBranchAddress("muppT", &entry[4]);
    chain.SetBranchAddress("mupy", &entry[5]);
    chain.SetBranchAddress("mumE", &entry[6]);
    chain.SetBranchAddress("mumpT", &entry[7]);
    chain.SetBranchAddress("mumy", &entry[8]);


    for (int i(0); i < chain.GetEntries(); ++i) {
        chain.GetEntry(i);

        for (int j(0); j < 9; ++j) {
            JpsiHist[j]->Fill(entry[j]);
        }
    }

    //Plot settings
    gStyle->SetOptStat(0);
    gStyle->SetLegendTextSize(0.03);
    TGaxis::SetMaxDigits(4);            //force exponential notation on axes

    //Make cuts here 
    TCanvas* cuts = new TCanvas("cuts", "cuts", 500, 500);

    //draw rapidity ditribution
    JpsiHist[2]->SetTitle("J/#psi Rapidity with Acceptance Cuts; y [GeV/c]; Events");
    JpsiHist[2]->SetMarkerStyle(21);
    JpsiHist[2]->SetMarkerColor(4);
    JpsiHist[2]->Draw("P");

    //preform cuts on acceptance
    TH1F* ycut1 = new TH1F("ycut1", "ycut1", numBins, -6.5, 6.5);
    ycut1->SetMarkerStyle(22);
    ycut1->SetMarkerColor(2);
    // Cuts in the form: 
    // mu+ -> y cut, pT cut
    // mu- -> y cut, pT cut
    //chain.Draw("Jpsiy>>ycut1", " \
        abs(mumy) < 1.2 && abs(mumpT) > 3.5 && \
        abs(mupy) < 1.2 && abs(muppT) > 3.5 || \
        1.2 < abs(mumy) && abs(mumy) < 2.1 && abs(mumpT) > (5.47 - 1.89 * abs(mumy)) && \
        1.2 < abs(mupy) && abs(mupy) < 2.1 && abs(muppT) > (5.47 - 1.89 * abs(mupy)) || \
        2.1 < abs(mumy) && abs(mumy) < 2.4 && abs(mumpT) > 1.5 && \
        2.1 < abs(mupy) && abs(mupy) < 2.4 && abs(muppT) > 1.5 \
        ", "P same");
    //go from 
    chain.Draw("Jpsiy>>ycut1", "abs(Jpsiy) < 0.5", "P same");
    //make selection in only positive muon
    //chain.Draw("Jpsiy>>ycut1", " \
        abs(mumy) < 1.2 && abs(mumpT) > 3.5 || \
        1.2 < abs(mumy) && abs(mumy) < 2.1 && abs(mumpT) > (5.47 - 1.89 * abs(mumy)) || \
        2.1 < abs(mumy) && abs(mumy) < 2.4 && abs(mumpT) > 1.5 \
        ", "P same");

    //make legend
    TLegend* leg = new TLegend(0.5, 17200, 7.5, 19675, "", "");
    leg->AddEntry(JpsiHist[2], "Complete y distribution");
    leg->AddEntry(ycut1, "Cuts from CMS");
    leg->Draw("same");

    TLatex TeX;
    TeX.DrawLatex(-7.0, 17700, Form("#scale[0.6]{#splitline{Complete y to CMS y }{%.0f/%.0f = %.3f}}", \
        ycut1->Integral(), JpsiHist[2]->Integral(), ycut1->Integral() / JpsiHist[2]->Integral()));

    //make plot for pT, y space
    TCanvas* pty = new TCanvas("pty", "pty", 500, 500);

    //make 2D histogram 
    TH2D* ptyhist = new TH2D("ptyhist", "J/#psi p_{T}, y Space with Acceptance Cuts; y; p_{T}; events", \
        numBins, -2.7, 2.7, numBins, 0.0, 30);

    chain.Draw("JpsipT:Jpsiy>>ptyhist", " \
        abs(mumy) < 1.2 && abs(mumpT) > 3.5 && \
        abs(mupy) < 1.2 && abs(muppT) > 3.5 || \
        1.2 < abs(mumy) && abs(mumy) < 2.1 && abs(mumpT) > (5.47 - 1.89 * abs(mumy)) && \
        1.2 < abs(mupy) && abs(mupy) < 2.1 && abs(muppT) > (5.47 - 1.89 * abs(mupy)) || \
        2.1 < abs(mumy) && abs(mumy) < 2.4 && abs(mumpT) > 1.5 && \
        2.1 < abs(mupy) && abs(mupy) < 2.4 && abs(muppT) > 1.5 \
        ", "colz");

    //chain.Draw("JpsipT:Jpsiy>>ptyhist", " \
        abs(mumy) < 1.2 && abs(mumpT) > 3.5 || \
        1.2 < abs(mumy) && abs(mumy) < 2.1 && abs(mumpT) > (5.47 - 1.89 * abs(mumy)) || \
        2.1 < abs(mumy) && abs(mumy) < 2.4 && abs(mumpT) > 1.5 \
        ", "colz");

}