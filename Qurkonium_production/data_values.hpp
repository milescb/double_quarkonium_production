// 
//  Miles Cochran-Branson, 8/2/21
//
//  Values below to be used in multinomial model 
//  for double production of quarkonium. 
//  Values come from PYTHIA and CMS. Calculation scrips can 
//  be found in pythia_simulations folder. 
//

class DataVals {

    private:

        //Calculated values from data

        //PARAMETERS OF MODEL

        //resaerched values here
        double ppCrossSection = 70; //units: mb
        double upsilonBR = 0.025; //2.5% branching ratio
        double jpsiBR = 0.06; //6% branching ratio
        //double upsilonBR = 1; //2.5% branching ratio
        //double jpsiBR = 1; //6% branching ratio
        double upsscaleFactor = 0.15; //go from one unit of rapidity to go to all, upsilon
        //double jpsiscaleFactor = 0.117; //go from one unit of rap to all rap, J/psi
        double jpsiscaleFactor = 0.121; //new values!! 
        double upsilonFdirect = 0.66; // (60 +/- 10)%
        double jpsiFdirect = 0.6; // (66 +/- 10)%

        //PROBABILITIES

        //cross sections from PYTHIA
        double upsilonCrossSection1 = 5.736e-04; //units: mb
        double jpsiCrossSection1 = 2.752e-2; //units: mb
        //modify quoted cross section appropriately 
        double upsilonCrossSection_final1 = upsilonBR * upsilonCrossSection1 / upsilonFdirect;
        double jpsiCrossSection_final1 = jpsiBR * jpsiCrossSection1 / jpsiFdirect;

        //report probability 
        double UpsilonProb1 = upsilonCrossSection_final1 / ppCrossSection;
        double JpsiProb1 = jpsiCrossSection_final1 / ppCrossSection;

        //cross sections from Feng et. al
        double upsilonCrossSection2 = 931e-09; //units: mb
        double jpsiCrossSection2 = 192.3e-06; //units: mb
        //modify quoted cross section appropriately 
        double upsilonCrossSection_final2 = upsilonCrossSection2 / (upsilonFdirect * upsscaleFactor);
        double jpsiCrossSection_final2 = jpsiCrossSection2 / (jpsiFdirect * jpsiscaleFactor);

        //report probability 
        double UpsilonProb2 = upsilonCrossSection_final2 / ppCrossSection;
        double JpsiProb2 = jpsiCrossSection_final2 / ppCrossSection;

        //cross sections from d'Enterria and Snigirev
        //scale from 5.5 TeV to 5.02 TeV; scale comes from curves from Feng et al paper
        double ups55to502scale = 909.6 / 976.6; 
        double jpsi55to502scale = 192.3 / 209.117;

        //scale from 5.02 TeV to 8.16 TeV (for PbPb run to pPb run)
        // double ups55to502scale = 1896 / 976.6; 
        // double jpsi55to502scale = 255 / 209.117;

        double jpsiCrossSection3 = 25e-03; //units: mb
        double upsilonCrossSection3 = 1.7e-03; //units: mb
        //modify quoted cross section appropriately 
        //is branching ratio used here? additionally: is this prompt or direct? 
        double upsilonCrossSection_final3 = upsilonCrossSection3 * ups55to502scale;
        double jpsiCrossSection_final3 = jpsiCrossSection3 * jpsi55to502scale;

        //report probability 
        double UpsilonProb3 = upsilonCrossSection_final3 / ppCrossSection;
        double JpsiProb3 = jpsiCrossSection_final3 / ppCrossSection;

        //ACCEPTANCE AND EFFICIENCY

        //acceptance calculated from PYTHIA
        double UpsilonAcceptance = 0.25;
        double JpsiAcceptance = 0.014;
        //efficiency from convolution
        double UpsilonEfficiency = 0.80;
        double JpsiEfficiency = 0.67;
        //Combined efficiency and acceptance from calculations
        double UpsilonA_E = 0.196918;
        double JpsiA_E = 0.0093773;

        //NUMBER OF EVENTS

        //at sqrt(s) = 5.02 TeV
        double luminosity = 1.7 * 1e+9; //units: inverse barn
        double PbPbcrosssection = 7.52191; //units: barn
        double NEvents = luminosity * PbPbcrosssection;

        //at sqrt(s) = 8.16 TeV
        double lum2 = 180 * 1e+9; //units: inverse barn
        double pPbcrosssection = 2.09; //units: barn
        //double NEvents = lum2 * pPbcrosssection;

    public:

        //At instance of class ask which probability to use 
        //Ask if A X E will be used
        int whichProb;
        bool useAcceptEff;
        DataVals() {
            cout << "Which probabilty? (1 = PYTHIA, 2 = Feng, 3 = d'Enterrria) ";
            cin >> whichProb;
            cout << "Use Acceptance + Efficiency? ";
            cin >> useAcceptEff;
        }

        //return probability of Upsilon
        double UpsProb() {
            if (whichProb == 1) return UpsilonProb1;
            else if (whichProb == 2) return UpsilonProb2;
            else if (whichProb == 3) return UpsilonProb3;
            else return 0.0;
        }

        //return probability of Jpsi
        double JpsiProb() {
            if (whichProb == 1) return JpsiProb1;
            else if (whichProb == 2) return JpsiProb2;
            else if (whichProb == 3) return JpsiProb3;
            else return 0.0;
        }

        //return Acceptance and efficiency 
        double UpsAE() {return UpsilonA_E;}
        double JpsiAE() {return JpsiA_E;}

        //return number of total events
        double NumEvents() {return NEvents;}

        //To access intermediate values for printing purposes
        double UpsProb1() {return UpsilonProb1;}
        double UpsProb2() {return UpsilonProb2;}
        double UpsProb3() {return UpsilonProb3;}
        double jpsiProb1() {return JpsiProb1;}
        double jpsiProb2() {return JpsiProb2;}
        double jpsiProb3() {return JpsiProb3;}
        double UpsCross1() {return upsilonCrossSection_final1;}
        double UpsCross2() {return upsilonCrossSection_final2;}
        double UpsCross3() {return upsilonCrossSection_final3;}
        double JpsiCross1() {return jpsiCrossSection_final1;}
        double JpsiCross2() {return jpsiCrossSection_final2;}
        double JpsiCross3() {return jpsiCrossSection_final3;}


};
////////////////////
//pPb stats 
////////////////////

// //cross sections from Feng et. al
// double upsilonCrossSection2 = 1896e-09; //units: mb
// double jpsiCrossSection2 = 255e-06; //units: mb
// //modify quoted cross section appropriately 
// double upsilonCrossSection_final2 = upsilonCrossSection2 / (upsilonFdirect * upsscaleFactor);
// double jpsiCrossSection_final2 = jpsiCrossSection2 / (jpsiFdirect * jpsiscaleFactor);

// //report probability 
// double UpsilonProb2 = upsilonCrossSection_final2 / ppCrossSection;
// double JpsiProb2 = jpsiCrossSection_final2 / ppCrossSection;

// //pPb stats
// double lum2 = 180 * 1e+9; //units: inverse barn
// double pPbcrosssection = 2.09; //units: barn
// double NEvents = lum2 * pPbcrosssection;