//
//  glauber.cpp
//  Glauber model simulation base code!
//
//  Created by Miles Cochran-Branson on 6/23/21.
//
#include <iostream>
#include <TNtuple.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom3.h>

using namespace std;

//Initial values for Pb nucleous
int A = 208;
//double R = 1.07 * pow(A,(1.0/3.0)); //differs here from actual
double a = 0.546;
double R = 6.62; //this is *actually* the value for Pb207... what's correct??!!

//define volumn density function
double rho(double r) {
    return 1 / (1 + exp((r - R) / a));
}

//Define probability density
double rho_dist(double r) {
    return ((r * r) / 25.0) / (1 + exp((r - R) / a));
}

//define centrality thang
double impact_param(double b) {
    return 2 * M_PI * b;
}

//input values from fit
double k = 0.746838;
double mu = 1.42448;

//define negative bionomial dist with gamm function
double nbd(double n) {
    return (tgamma(n + k) / (tgamma(n + 1) * tgamma(k))) * (pow((mu / k), n) / (pow((mu / k) + 1, (n + k))));
}

//define number of collisions
//int numNuclei = 1e6;

void glauber(int numNuclei = 10) {
    
    //***INITIALIZE OUTPUT STATEMENTS***//
    
    //create counter for print statements
    int counter(0);
    
    //Create file to output results!
    TFile ofile("glauberfinal1.root", "RECREATE");
    
    //Create NTuple
    TNtuple glauberfinal1("glauberfinal1", "Ntuple Data", "ImpactParameter:NumCollisions:NumParticipating:SumEt");
    
    //create function to pass through rng
    TF1* rhoDist = new TF1("rhoDist", "rho_dist(x)", 0.0, 10.0);
    
    //Define function for random impact param
    TF1 *impact = new TF1("impact", "impact_param(x)",0.0, 18.0);
    
    //Create TF1 from which to sample random numbers
    TF1* nbdist = new TF1("nbdist", "nbd(x)", 0.0, 20.0);
    
    //***BEGIN BULK OF COMPUTATION***//
    
    for (int iNucleon(0); iNucleon < numNuclei; ++iNucleon) {
        
        //****BEGIN POPULATING NUCLEI****//
        
        //Initialize random number generator
        TRandom3 rnd3(0);
            
        //Create array to hold nucleon possition
        double x1[A], y1[A], z1[A];
        
        //Create array to hold second nucleon possition
        double x2[A], y2[A], z2[A];
        
        //Loop over all nucleons for first nucleus
        for (int iNucleons = 0; iNucleons < A; ++iNucleons) {
        
            //Create random dist. in spherical coords.
            double radius1 = rhoDist -> GetRandom();
            double radius2 = rhoDist -> GetRandom();
            double theta1 = acos(rnd3.Uniform(-1.0, 1.0));
            double theta2 = acos(rnd3.Uniform(-1.0, 1.0));
            double phi1 = rnd3.Uniform(0.0, 2 * M_PI);
            double phi2 = rnd3.Uniform(0.0, 2 * M_PI);
        
            //transform to cartesian
            x1[iNucleons] = radius1 * cos(phi1) * sin(theta1);
            y1[iNucleons] = radius1 * sin(phi1) * sin(theta1);
            z1[iNucleons] = radius1 * cos(theta1);
            x2[iNucleons] = radius2 * cos(phi2) * sin(theta2);
            y2[iNucleons] = radius2 * sin(phi2) * sin(theta2);
            z2[iNucleons] = radius2 * cos(theta2);
        }
        
        //****END NUCLEI POPULATION//
        
        //****BEGIN SET RANDOM OFFSET PARAMETER***//
        
        //Set Impact parameter
        double impactParamOffset = impact->GetRandom();
        double impactParam[A];
        for (int i = 0; i < A; i++) { impactParam[i] = impactParamOffset; } //make impact param into array
        
        //offset one nucleus on the x-axis by impact parameter b
        for (int i = 0; i < A; i++) { x2[i] = x2[i] - impactParam[i]; }
        
        //****END OFFSET PARAMETER****//
        
        //****BEGIN CALCULATION OF COLLISIONS/PARTICIPANTS****//
        
        //Set cross-sectional area, sigma
        double sigma(70.0); //units: mb  70mb sigma +/- 5mb for sqrt(s) = 5.02TeV collisions
        sigma = sigma / 10.0; // convert to fm^2
        
        //Create thing for nuber of collisions
        int numCollisions(0);
        
        //Create sets to store indices of par
        set<int>nucleus1Index;
        set<int>nucleus2Index;
        
        //create value for convolution
        double SumEt(0);
        
        //Now check who is colliding / not — find the distance
        for (int iNucleon1(0); iNucleon1 < A; ++iNucleon1) {
            
            bool isParticipating(0);
            
            for (int iNucleon2(0); iNucleon2 < A; ++iNucleon2) {
                
                double x_dis, y_dis, distance_sq;
                
                x_dis = x1[iNucleon1] - x2[iNucleon2];
                y_dis = y1[iNucleon1] - y2[iNucleon2];
                
                distance_sq = pow(x_dis, 2) + pow(y_dis, 2);
                
                if (distance_sq < (sigma / M_PI)) {
                    
                    numCollisions++;
                    isParticipating = 1;
                    nucleus2Index.insert(iNucleon2);
                    
                    //for convolution sample from negative binonmial
                    SumEt += nbdist->GetRandom();
                    
                }
                
            if (isParticipating == 1) {
                
                //numParticipating1++;
                nucleus1Index.insert(iNucleon1);
                
                }
            }
        }
        
        int numParticipating = nucleus1Index.size() + nucleus2Index.size();
        
        //****END CALCULATIONS****//
        
        
        //****PRINT RESULTS!****//
        
        //Output value
        
        if (iNucleon % 100000 == 0) {
            cout << "Impact Parameter: " << impactParamOffset << endl;
            cout << "Number of Collisions: " << numCollisions << endl;
            cout << "Number of Participants: " << numParticipating << endl;
            cout << "Counter value: " << counter << endl;
            cout << "NEW NUCLEI" << endl;
            ++counter;
        }
        
        //Write data to our NTuple
        if (numCollisions != 0) {
            glauberfinal1.Fill(impactParamOffset, numCollisions, numParticipating, SumEt / 1000);
        }
    }
    
    glauberfinal1.Write();
    ofile.Close();
    
    return;
}

int main() {
   glauber();
   return 0;
}
