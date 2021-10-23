// Author: Miles Cochran-Branson, 7/26/21
//
// cuts for Upsilon data. Requires .root file created from PYTHIA simulation
// makes kinematic cuts to determine acceptance and efficiency 

void upsilonCuts() {

    //extract info from root files
    TChain chain("upsDat");
    chain.Add("upsilonData_final1*.root");

    int numBins = 100;

    //make histograms here 
    TH1F* upsilonHist[9];
    int k(0);
    upsilonHist[k++] = new TH1F("upsE", "#varUpsilon Energy; Energy [GeV]; Events", numBins, -0.5, 1999.5);
    upsilonHist[k++] = new TH1F("upspT", "#varUpsilon p_{T}; p_{T} [GeV/c]; Events", numBins, -0.5, 30);
    upsilonHist[k++] = new TH1F("upsRap", "#varUpsilon Rapidity; y; Events", numBins, -6.5, 6.5);
    upsilonHist[k++] = new TH1F("muonpE", "#mu^{+} Energy; Energy [GeV]; Events", numBins, -0.5, 1999.5);
    upsilonHist[k++] = new TH1F("muonppT", "#mu^{+} p_{T}; p_{T} [GeV/c]; Events", numBins, -0.5, 30);
    upsilonHist[k++] = new TH1F("muonpy", "#mu^{+} Rapidity; y; Events", numBins, -6.5, 6.5);
    upsilonHist[k++] = new TH1F("muonmE", "#mu^{-} Energy; Energy [GeV]; Events", numBins, -0.5, 1999.5);
    upsilonHist[k++] = new TH1F("muonmpT", "#mu^{-} p_{T}; p_{T} [GeV/c]; Events", numBins, -0.5, 30);
    upsilonHist[k++] = new TH1F("muonmy", "#mu^{-} Rapidity; y; Events", numBins, -6.5, 6.5);
    
    //Set Branch Addresses 
    float entry[9];
    chain.SetBranchAddress("upsE", &entry[0]);
    chain.SetBranchAddress("upspT", &entry[1]);
    chain.SetBranchAddress("upsy", &entry[2]);
    chain.SetBranchAddress("mupE", &entry[3]);
    chain.SetBranchAddress("muppT", &entry[4]);
    chain.SetBranchAddress("mupy", &entry[5]);
    chain.SetBranchAddress("mumE", &entry[6]);
    chain.SetBranchAddress("mumpT", &entry[7]);
    chain.SetBranchAddress("mumy", &entry[8]);


    for (int i(0); i < chain.GetEntries(); ++i) {
        chain.GetEntry(i);

        for (int j(0); j < 9; ++j) {
            upsilonHist[j]->Fill(entry[j]);
        }
    }

    //Plot settings
    gStyle->SetOptStat(0);
    gStyle->SetLegendTextSize(0.03);
    TGaxis::SetMaxDigits(4);            //force exponential notation on axes

    //Make cuts here 
    TCanvas* cuts = new TCanvas("cuts", "cuts", 500, 500);

    //draw rapidity ditribution
    upsilonHist[2]->SetTitle("#varUpsilon Rapidity with Acceptance Cuts; y [GeV/c]; Events");
    upsilonHist[2]->SetMarkerStyle(21);
    upsilonHist[2]->SetMarkerColor(4);
    upsilonHist[2]->Draw("P");

    //preform cuts on efficieny and acceptance
    TH1F* ycut1 = new TH1F("ycut1", "ycut1", numBins, -6.5, 6.5);
    TH1F* ycut2 = new TH1F("ycut2", "ycut2", numBins, -6.5, 6.5);
    ycut1->SetMarkerStyle(22);
    ycut1->SetMarkerColor(2);
    ycut2->SetMarkerStyle(23);
    ycut2->SetMarkerColor(6);
    chain.Draw("upsy>>ycut1", "abs(mupy) < 2.4 && abs(muppT) > 3.5 && \
        abs(mumy) < 2.4 && abs(mumpT) > 3.5", "P same");
    //chain.Draw("upsy>>ycut2", "abs(mupy) < 0.5 && \
        abs(mumy) < 0.5", "P same");
    //go from one unit of rap to all rap
    chain.Draw("upsy>>ycut2", "abs(upsy) < 0.5", "P same");
    //make legend
    TLegend* leg = new TLegend(0.7, 16750, 6.5, 20500, "", "");
    leg->AddEntry(upsilonHist[2], "Complete y distribution");
    leg->AddEntry(ycut1, "|y^{#mu^{+}#mu^{-}}| < 2.4, |p_{T}^{#mu^{+}#mu^{-}}| > 3.5");
    leg->AddEntry(ycut2, "|y^{#mu^{+}#mu^{-}}| < 0.5, |p_{T}^{#mu^{+}#mu^{-}}| > 3.5");
    leg->Draw("same");

    TLatex TeX;
    TeX.DrawLatex(-6.2, 19000, Form("#scale[0.5]{#splitline{Complete y to CMS y }{%.0f/%.0f = %.2f}}", \
        ycut1->Integral(), upsilonHist[2]->Integral(), 1/(upsilonHist[2]->Integral() / ycut1->Integral())));
    TeX.DrawLatex(-6.2, 17000, Form("#scale[0.5]{#splitline{Complete y to |y| < 0.5}{%.0f/%.0f = %.2f}}", \
        ycut2->Integral(), upsilonHist[2]->Integral(), 1/(upsilonHist[2]->Integral() / ycut2->Integral())));
    TeX.DrawLatex(-6.2, 15000, Form("#scale[0.5]{#splitline{|y| < 2.4 to |y| < 0.5}{%.0f/%.0f = %.2f}}", \
        ycut2->Integral(), ycut1->Integral(), 1/(ycut1->Integral() / ycut2->Integral())));

    //plot for pT vs y space
    TCanvas* pty = new TCanvas("pty", "pty", 500, 500);

    //make 2D histogram 
    TH2D* ptyhist = new TH2D("ptyhist", "#varUpsilon p_{T}, y Space with Acceptance Cuts; y; p_{T}; events", \
        numBins, -2.5, 2.5, numBins, 0.0, 30);

    //draw with cuts
    chain.Draw("upspT:upsy>>ptyhist", "abs(mupy) < 2.4 && abs(muppT) > 3.5 && \
        abs(mumy) < 2.4 && abs(mumpT) > 3.5", "colz");
}