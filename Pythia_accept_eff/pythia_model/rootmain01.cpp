//
//  Simulation of Upsilon production to determine basic kinematic variables for cuts
//  Must be run with the Makefile: requires flags for pythia and root
//
//  Created by Miles Cochran-Branson on 7/20/21.
//
// include pythia packages as well as root packages
#include <iostream>
#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TNtuple.h"

using namespace Pythia8;

//define number of events
int numEvents = 1e6;

void UpsilonProduction() {
    
    //initialize Pythia
    Pythia pythia;

    //settings for decays
    pythia.readString("PhaseSpace:pTHatMin = 1."); //transverse momentum 
    pythia.readString("Beams:eCM = 5020"); //sqrt(s) = 5.02TeV
    //pythia.readString("Bottomonium:all = on"); //turn on production of Bottomonium
    pythia.readString("Bottomonium::gg2bbbar(3S1)[3S1(1)]g = {on, on, on}");
    pythia.readString("Bottomonium::gg2bbbar(3S1)[3S1(1)]gm = {on, on, on}");
    pythia.readString("Bottomonium::gg2bbbar(3S1)[3S1(8)]g = {on, on, on}");
    pythia.readString("Bottomonium::qg2bbbar(3S1)[3S1(8)]q = {on, on, on}");
    pythia.readString("Bottomonium::qqbar2bbbar(3S1)[3S1(8)]g = {on, on, on}");
    pythia.readString("Bottomonium::gg2bbbar(3S1)[1S0(8)]g = {on, on, on}");
    pythia.readString("Bottomonium::qg2bbbar(3S1)[1S0(8)]q = {on, on, on}");
    pythia.readString("Bottomonium::qqbar2bbbar(3S1)[1S0(8)]g = {on, on, on}");

    //requirments for upsilon decays
    //allow Upsilon (553) to only decay to mu- (13), mu+ (-13) 
    pythia.readString("553:onMode = off");
    pythia.readString("553:onIfMatch = 13 -13");

    pythia.init();
    
    //create root output file
    TFile* ofile = new TFile("upsilonData.root", "RECREATE");
    
    //Ntuple to store data. Alias: ups = Upsilon, mup = mu-minus, mum = mu-plus;
    //E = energy, pT = transverse p, y = rapidity, phi = angle phi 
    TNtuple* upsDat = new TNtuple("upsDat", "data", "upsE:upspT:upsy:upsphi:upsm:mupE:muppT:mupy:mupphi:mumE:mumpT:mumy:mumphi");

    float nTupleArray[13];
    
    for (int iEvent(0); iEvent < numEvents; ++iEvent) {
        if (!pythia.next()) continue;
        
        for (int i(0); i < pythia.event.size(); ++i) {
            if (pythia.event[i].id() == 553) {

                //Locate and identify bottom location of Upsilon
                Particle& upsilonParticleArb = pythia.event[i];
                int iBotUpsilon = upsilonParticleArb.iBotCopy();
                Particle& upsilonParticle = pythia.event[iBotUpsilon];

                //get 4-momentum of Upsilon
                Vec4 ups4Mom = upsilonParticle.p();

                //locate and identify daughter particles
                int iDaughter1 = upsilonParticle.daughter1();
                int iDaughter2 = upsilonParticle.daughter2();
                //identify which is mu- / mu+
                int iMuPlus;
                int iMuMinus;
                if (pythia.event[iDaughter1].id() == 13) {
                    iMuMinus = iDaughter1;
                    iMuPlus = iDaughter2;
                }
                else if (pythia.event[iDaughter1].id() == -13) {
                    iMuMinus = iDaughter2;
                    iMuPlus = iDaughter1;
                }
                else {
                    cout << "Error: Decay is not dimuonic" << endl;
                    break;
                }
                Particle& MuMinus = pythia.event[iMuMinus];
                Particle& MuPlus = pythia.event[iMuPlus];
            
                //define 4-momentum of products
                Vec4 MuMinus4Mom = MuMinus.p();
                Vec4 MuPlus4Mom = MuPlus.p();

                //update array with events
                int k = 0;
                nTupleArray[k++] = ups4Mom.e();
                nTupleArray[k++] = ups4Mom.pT();
                nTupleArray[k++] = ups4Mom.rap();
                nTupleArray[k++] = ups4Mom.phi();
                nTupleArray[k++] = ups4Mom.mCalc();
                nTupleArray[k++] = MuPlus4Mom.e();
                nTupleArray[k++] = MuPlus4Mom.pT();
                nTupleArray[k++] = MuPlus4Mom.rap();
                nTupleArray[k++] = MuPlus4Mom.phi();
                nTupleArray[k++] = MuMinus4Mom.e();
                nTupleArray[k++] = MuMinus4Mom.pT();
                nTupleArray[k++] = MuMinus4Mom.rap();
                nTupleArray[k++] = MuMinus4Mom.phi();

                int percent = 100000;
                if (i % percent == 0) {

                    cout << Form("%f%% complete", percent * 100 / numEvents ) << endl;
                    cout << " " << endl;

                    //print statements for Upsilon
                    cout << "Upsilon found at location: " << iBotUpsilon << endl;
                    cout << "Upsilon 4-momentum " << ups4Mom << endl;
                
                    //print statements for Dimuon
                    cout << "Mu- found at location: " << MuMinus.index() << endl;
                    cout << "Mu- 4-momentum: " << MuMinus4Mom << endl;
    
                    cout << "Mu+ found at location: " << MuPlus.index() << endl;
                    cout << "Mu+ 4-momentum: " << MuPlus4Mom << endl;
                }

                //end loop after finding Upsilon and daughters
                break;
            }
            
        }

        //fill NTuple
        upsDat->Fill(nTupleArray);
        
    }
    
    //write then close file
    upsDat->Write();
    ofile->Close();

    pythia.stat();
    
    return;
}

int main() {
    UpsilonProduction();
    return 0;
}
