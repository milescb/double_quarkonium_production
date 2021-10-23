//
// Miles Cochran-Branson, 8/2/21
//
// Calculate efficiency from CMS data. 
// Requires:
// .root file with data from CMS
// .root file from PYTHIA
//  Acceptance cuts! as root file from CMS includes acceptance

//add histograms of pseudo-experiments

#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TNtuple.h>

using namespace std;

int numSimulations = 1000;

void efficiency_final() {

    //gStyle->SetOptFit(111);
    //gStyle->SetLineWidth(4);

    //read .root file from CMS
    TFile* HS_eff = new TFile("HS_eff.root");
    TH2D* mueff_initial = (TH2D*)HS_eff->Get("muAccEffEtaPt");

    TFile* ofile = new TFile("efficiencyAcceptanceData.root", "RECREATE");
    TNtuple* effAcceptData = new TNtuple("effAcceptData", "data", "upsAE:jpsiAE");

    //TCanvas* c = new TCanvas("c", "c", 500, 500);
    //mueff_initial->Draw("colz");

    int xBins = 2 * mueff_initial->GetNbinsX();
    int yBins = 5 * mueff_initial->GetNbinsY(); //multiply by 5 to scale to pT = 30; (30 / 6)

    //re-create histogram s.t we have data until 30. 
    TH2D* mueff = new TH2D("mueff", "#mu Efficiency and Acceptance Plot from CMS data; y; p_{T}", \
        xBins, -2.7, 2.7, yBins, 0.0, 30);

    //set bin contents of mueff
    double AvgLastContent; 
    int numAvgDown(5);
    // add reflect over x = 0 axis as well
    for (int i(0); i < xBins; ++i) {
        for (int j(0); j < yBins; ++j) {
            if (i <= (xBins / 2)) { //make mirror image of graph to include all rapidity space
               if (j > mueff_initial->GetNbinsY()) { //condition for pT > 6
                    AvgLastContent = 0;
                    for (int k(0); k < numAvgDown; ++k) {
                        AvgLastContent += mueff_initial->GetBinContent((xBins/2) - i + 1, mueff_initial->GetNbinsY() - k + 1);
                    }
                    AvgLastContent = AvgLastContent / numAvgDown;
                    mueff->SetBinContent(i, j, AvgLastContent);
                }
                else {
                    mueff->SetBinContent(i, j, mueff_initial->GetBinContent((xBins/2) - i + 1, j));
                } 
            }
            else {
                if (j > mueff_initial->GetNbinsY()) { //condition for pT > 6
                    AvgLastContent = 0;
                    for (int k(0); k < numAvgDown; ++k) {
                        AvgLastContent += mueff_initial->GetBinContent(i - (xBins/2), mueff_initial->GetNbinsY() - k + 1);
                    }
                    AvgLastContent = AvgLastContent / numAvgDown;
                    mueff->SetBinContent(i, j, AvgLastContent);
                }
                else {
                    mueff->SetBinContent(i, j, mueff_initial->GetBinContent(i - (xBins / 2), j));
                }
            }
        }
    }

    //draw CMS plot
    //TCanvas* c = new TCanvas("c", "c1", 500, 500);

    //mueff->Draw("colz");

    //****Read pythia data and re-create 2D histograms of data with cuts****//

    //****Upsilon data****//

    //extract info from root files
    TChain upschain("upsDat");
    upschain.Add("upsilonData_final1*.root");
    
    //Set Branch Addresses 
    float entryups[9];
    upschain.SetBranchAddress("upsE", &entryups[0]);
    upschain.SetBranchAddress("upspT", &entryups[1]);
    upschain.SetBranchAddress("upsy", &entryups[2]);
    upschain.SetBranchAddress("mupE", &entryups[3]);
    upschain.SetBranchAddress("muppT", &entryups[4]);
    upschain.SetBranchAddress("mupy", &entryups[5]);
    upschain.SetBranchAddress("mumE", &entryups[6]);
    upschain.SetBranchAddress("mumpT", &entryups[7]);
    upschain.SetBranchAddress("mumy", &entryups[8]);

    //****J/psi data****//

    //extract info from root files
    TChain jpsichain("JpsiDat");
    jpsichain.Add("JpsiData_final1*.root");
    
    //Set Branch Addresses 
    float entryjpsi[9];
    jpsichain.SetBranchAddress("JpsiE", &entryjpsi[0]);
    jpsichain.SetBranchAddress("JpsipT", &entryjpsi[1]);
    jpsichain.SetBranchAddress("Jpsiy", &entryjpsi[2]);
    jpsichain.SetBranchAddress("mupE", &entryjpsi[3]);
    jpsichain.SetBranchAddress("muppT", &entryjpsi[4]);
    jpsichain.SetBranchAddress("mupy", &entryjpsi[5]);
    jpsichain.SetBranchAddress("mumE", &entryjpsi[6]);
    jpsichain.SetBranchAddress("mumpT", &entryjpsi[7]);
    jpsichain.SetBranchAddress("mumy", &entryjpsi[8]);

    //initialize rng
    TRandom3 rnd3(0);
    double rand1, rand2;

    //now we perform many Monte-Carlo simulations 

    //upsilon stats
    int upsMuPlusBin, upsMuMinusBin;
    double upsMuPlusVal, upsMuMinusVal; 
    int totUpsilon;
    int acceptUpsilon;
    int reconstructedUpsilon;

     //jpsi stats here
    int jpsiMuPlusBin, jpsiMuMinusBin;
    double jpsiMuPlusVal, jpsiMuMinusVal; 
    int totJpsi;
    int acceptJpsi;
    int reconstructedJpsi;

    double mumy;
    double mumpT;
    double mupy;
    double muppT;

    double upsilonAcceptEfficiency = 0.0;
    double upsilonAcceptance = 0.0;

    double jpsiAcceptEfficiency = 0.0;
    double jpsiAcceptance = 0.0;

    for (int inumSim(0); inumSim < numSimulations; ++inumSim) {

        //print status in 1% increments :
        if (inumSim % 100 == 0) {
            cout << Form("%.1f%% complete", static_cast<double>(inumSim) * 100.0 / static_cast<double>(numSimulations)) << endl;
        }

        totUpsilon = 0;
        reconstructedUpsilon = 0;
        acceptUpsilon = 0;

        totJpsi = 0;
        reconstructedJpsi = 0;
        acceptJpsi = 0;

        //loop over all data in the both ntuples
        //check if both muons are in the effienciy defined by CMS data
        for (int i(0); i < jpsichain.GetEntries(); ++i) {

            rand1 = rnd3.Rndm();
            rand2 = rnd3.Rndm();

            //upsilon here
            upschain.GetEntry(i);

            //increment total Upsilon
            ++totUpsilon;

            //make acceptance cuts for Upsilon here
            if ((abs(entryups[5]) < 2.4) && (abs(entryups[4]) > 3.5) && (abs(entryups[8]) < 2.4) && (abs(entryups[7]) > 3.5)) {

                ++acceptUpsilon;

                upsMuPlusBin = mueff->FindBin(entryups[5], entryups[4]);
                upsMuPlusVal = mueff->GetBinContent(upsMuPlusBin);

                upsMuMinusBin = mueff->FindBin(entryups[8], entryups[7]);
                upsMuMinusVal = mueff->GetBinContent(upsMuMinusBin);

                if (upsMuPlusVal > rand1 && upsMuMinusVal > rand2) ++reconstructedUpsilon;
            }

            //jpsi data
            jpsichain.GetEntry(i);

            mumy = entryjpsi[8];
            mumpT = entryjpsi[7];
            mupy = entryjpsi[5];
            muppT = entryjpsi[4];

            //increment total J/psi
            ++totJpsi;

            //make accpetance cuts for J/psi here
            if ( ((abs(mumy) < 1.2 && abs(mumpT) > 3.5) && \
                    ((abs(mupy) < 1.2 && abs(muppT) > 3.5) || (1.2 < abs(mupy) && abs(mupy) < 2.1 && abs(muppT) > (5.47 - 1.89 * abs(mupy))) || \
                    (2.1 < abs(mupy) && abs(mupy) < 2.4 && abs(muppT) > 1.5))) || \
                ( (1.2 < abs(mumy) && abs(mumy) < 2.1 && abs(mumpT) > (5.47 - 1.89 * abs(mumy))) && \
                    ((abs(mupy) < 1.2 && abs(muppT) > 3.5) || (1.2 < abs(mupy) && abs(mupy) < 2.1 && abs(muppT) > (5.47 - 1.89 * abs(mupy))) || \
                    (2.1 < abs(mupy) && abs(mupy) < 2.4 && abs(muppT) > 1.5))) || \
                ( (2.1 < abs(mumy) && abs(mumy) < 2.4 && abs(mumpT) > 1.5) && \
                    ((abs(mupy) < 1.2 && abs(muppT) > 3.5) || (1.2 < abs(mupy) && abs(mupy) < 2.1 && abs(muppT) > (5.47 - 1.89 * abs(mupy))) || \
                    (2.1 < abs(mupy) && abs(mupy) < 2.4 && abs(muppT) > 1.5))) ) {
                
                ++acceptJpsi;

                jpsiMuPlusBin = mueff->FindBin(entryjpsi[5], entryjpsi[4]);
                jpsiMuPlusVal = mueff->GetBinContent(jpsiMuPlusBin); 

                jpsiMuMinusBin = mueff->FindBin(entryjpsi[8], entryjpsi[7]);
                jpsiMuMinusVal = mueff->GetBinContent(jpsiMuMinusBin);

                if (jpsiMuPlusVal > rand1 && jpsiMuMinusVal > rand2) ++reconstructedJpsi;  
            }
        }

        upsilonAcceptEfficiency += static_cast<double>(reconstructedUpsilon) / static_cast<double>(totUpsilon);
        //ups->Fill(static_cast<double>(reconstructedUpsilon) / static_cast<double>(totUpsilon));
        upsilonAcceptance += static_cast<double>(acceptUpsilon) / static_cast<double>(totUpsilon);

        jpsiAcceptEfficiency += static_cast<double>(reconstructedJpsi) / static_cast<double>(totJpsi);
        //jpsi->Fill(static_cast<double>(reconstructedJpsi) / static_cast<double>(totJpsi));
        jpsiAcceptance += static_cast<double>(acceptJpsi) / static_cast<double>(totJpsi);

        effAcceptData->Fill(static_cast<double>(reconstructedUpsilon) / static_cast<double>(totUpsilon), \
            static_cast<double>(reconstructedJpsi) / static_cast<double>(totJpsi));
    }

    //report values of efficiency and acceptance here:

    double upsilonFinalAcceptEff = upsilonAcceptEfficiency / static_cast<double>(numSimulations);
    double upsilonAccept = upsilonAcceptance / static_cast<double>(numSimulations);
    double upsilonEfficiency = upsilonFinalAcceptEff / upsilonAccept;

    double jpsiFinalAcceptEff = jpsiAcceptEfficiency / static_cast<double>(numSimulations);
    double jpsiAccept = jpsiAcceptance / static_cast<double>(numSimulations);
    double jpsiEfficiency = jpsiFinalAcceptEff / jpsiAccept;

    cout << " " << endl;
    cout << "Upsilon acceptance = " << upsilonAccept << endl;
    cout << "Upsilon efficiency = " << upsilonEfficiency << endl;
    cout << "Upsilon efficiency + acceptance = " << upsilonFinalAcceptEff << endl;

    cout << " " << endl;
    cout << "J/psi acceptance = " << jpsiAccept << endl;
    cout << "J/psi efficiency = " << jpsiEfficiency << endl;
    cout << "J/psi efficiency + acceptance = " << jpsiFinalAcceptEff << endl;

    effAcceptData->Write();
    ofile->Close();

    return;

}

int main() {
    efficiency_final();
}