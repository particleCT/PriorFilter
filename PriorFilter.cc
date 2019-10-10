#include <iostream>
#include <fstream>
#include <math.h>       /* sqrt */
#include <cmath>        // std::abs
#include <stdlib.h>  
#include "TVector3.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph2D.h"
#include "TProfile2D.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;
struct Particle{
  Float_t x0,y0,z0,px0,py0,pz0;
  Float_t x1,y1,z1,px1,py1,pz1;
  Float_t wepl;
  Int_t   Id;
  Float_t WET_prob, Y_prob, Z_prob;
  Float_t angle;
};

std::string* part_name;
double EstimateEstop(Particle*, TH3D*, TMatrixD&, TMatrixD&, TMatrixD&, TMatrixD&);
double EnergyStraggling(TMatrixD, double);


double simps(TMatrixD , TMatrixD );
double E2beta(double );
vector<double> Energy;
vector<double> dEdXBins;

double findWET(double, double);
double findEnergy(double, double);
void simps_mult(TMatrixD, TMatrixD, double, double &, double &, double &);
double Gauss(double, double, double);
double Gauss2D(TMatrixD, TMatrixD, TMatrixD);
TMatrixD Mult(TMatrixD, TMatrixD);
TMatrixD Mult(TMatrixD, double );

int NStep = 500;
double det_Xmin = -16.43; // cm
double det_Xmax =  16.43; // cm
double det_Zmin = -8.76/2; // cm
double det_Zmax =  8.76/2; // cm
double det_Ymin = -35.2/2; // cm
double det_Ymax =  35.2/2; // cm

double reco_Xmin = -24.00/2.; // cm
double reco_Xmax =  24.00/2.; // cm
double reco_Ymin = -24.00/2.; // cm
double reco_Ymax =  24.00/2.; // cm
double reco_Zmin = -9.00/2.; // cm
double reco_Zmax =  9.00/2.; // cm


double midX = (det_Xmax+det_Xmin)/2.; // cm
double midY = (det_Ymax+det_Ymin)/2.; // cm
double midZ = (det_Zmax+det_Zmin)/2.; // cm

int main(int argc, char** argv){

  if(argc<=2) {
    cout<<"Please input the following arguments : Prior-File  Projection-File "<<endl;
    return 0;
  }

  Particle Point;
  //--------------------------------------
  // Load Phantom data
  char* phantomFileName = Form("%s",argv[1]);
  TFile* phantomFile = new TFile(phantomFileName,"update");
  char* histName = Form("RSP");
  TH3D* RSPMap = (TH3D*)phantomFile->Get(histName);

  //--------------------------------------
  // Load Particle data  
  char * phaseFilename = argv[2];
  TFile* f = new TFile(phaseFilename,"update");
  TTree* t = (TTree*)f->Get("phase");

  // get the data
  t->SetBranchAddress("x0",&Point.x0);
  t->SetBranchAddress("y0",&Point.y0);
  t->SetBranchAddress("z0",&Point.z0);

  t->SetBranchAddress("x1",&Point.x1);
  t->SetBranchAddress("y1",&Point.y1);
  t->SetBranchAddress("z1",&Point.z1);

  t->SetBranchAddress("px0",&Point.px0);
  t->SetBranchAddress("py0",&Point.py0);
  t->SetBranchAddress("pz0",&Point.pz0);

  t->SetBranchAddress("px1",&Point.px1);
  t->SetBranchAddress("py1",&Point.py1);
  t->SetBranchAddress("pz1",&Point.pz1);

  t->SetBranchAddress("angle",&Point.angle);
  t->SetBranchAddress("wepl",&Point.wepl);
  t->GetEntry(0);
  TBranch* bpt_WET  = t->Branch("WET_prob",&Point.WET_prob,"WET_prob/F");
  TBranch* bpt_Y    = t->Branch("Y_prob",&Point.Y_prob,"Y_prob/F");
  TBranch* bpt_Z    = t->Branch("Z_prob",&Point.Z_prob,"Z_prob/F");

  std::string line;
  std::ifstream SPWater ("dEdX/Water_Geant4.dat");

  double data[3];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<3;i++) ss >> data[i];
    Energy.push_back(data[0]);
    dEdXBins.push_back(data[1]);
  }
   
  // Do the filters
  TMatrixD R0(2,2);
  TMatrixD Y0(2,1);
  TMatrixD Z0(2,1);
  TMatrixD Y1(2,1);
  TMatrixD Z1(2,1);
  TMatrixD Sigma(2,2);
  double sigma_t , sigma_theta , sigma_t_theta,sigma_E, var_E, var_WET, sigma_WET;

  int NEntries = t->GetEntries();
  for(int i=0;i< NEntries;i++){

    t->GetEntry(i);

    if(i%10000==0) cout<<i<<endl;

    TMatrixD vector_X(NStep,1); // position in cm
    TMatrixD vector_L(NStep,1); // radiation length in cm
    TMatrixD vector_pv(NStep,1); // inverse energy
    TMatrixD vector_E(NStep,1); // energy     
    
    TVector3 p0(Point.px0,Point.py0,Point.pz0);
    TVector3 p1(Point.px1,Point.py1,Point.pz1);
    p0.SetMag(1);
    p1.SetMag(1);
    
    // Calculated the position expected value all in cm
    R0(0,0)  = 1        ; R0(0,1)  = (Point.x1 - Point.x0)/10;
    R0(1,0)  = 0        ; R0(1,1)  = 1;      
    Y0(0,0)  = Point.y0/10 ; Y0(1,0)   = p0.y() ;
    Z0(0,0)  = Point.z0/10 ; Z0(1,0)   = p0.z() ; 
    Y1(0,0)  = Point.y1/10 ; Y1(1,0)   = p1.y() ; 
    Z1(0,0)  = Point.z1/10 ; Z1(1,0)   = p1.z() ;

    double WEPL_mes  = Point.wepl/10;
    double Estop_mes = findEnergy(200,WEPL_mes);
    double Estop_est = EstimateEstop(&Point, RSPMap, vector_L, vector_pv, vector_E, vector_X);
    double WEPL_est  = findWET(200, Estop_est);

    // Calculate the standard deviation around the expected position

    double Epsilon   = pow(13.6,2);      
    double E0s0      = Epsilon*pow(1+0.038* TMath::Log(simps(vector_X, Mult(vector_L, R0(0,1) ) ) ),2);
    TMatrixD y1      = Mult( vector_pv, vector_L);
    simps_mult(vector_X, y1, vector_X(NStep-1,0), sigma_t, sigma_t_theta, sigma_theta);

    sigma_t       = E0s0*sigma_t;
    sigma_theta   = E0s0*sigma_theta;
    sigma_t_theta = E0s0*sigma_t_theta;

    // sttraging variance of energy -- This neglects scattering!
    var_E         = EnergyStraggling(vector_E, Estop_mes); 
    sigma_E       = TMath::Sqrt(var_E);

    // straggling variance of the WET
    int it_Estop  = lower_bound(Energy.begin(), Energy.end(), Estop_mes)-Energy.begin();
    double SW_estop  = dEdXBins[it_Estop];
    var_WET      = var_E/pow(SW_estop,2);
    sigma_WET    = TMath::Sqrt(var_WET);

    // Scattering matrix
    Sigma(0,0)    = sigma_t ;      Sigma(0,1) = sigma_t_theta;
    Sigma(1,0)    = sigma_t_theta; Sigma(1,1) = sigma_theta;

    // Minimal limit of 0.5 cm on the WET uncertainty -- this is detector
    sigma_WET< 0.5 ? sigma_WET =0.5: sigma_WET = sigma_WET;
    
    // Minimal limit of 0.5 cm on the scattering uncertainty
    sigma_t< 0.5 ? sigma_t =0.5: sigma_t = sigma_t;

    double Y_pred     = (R0*Y0)(0,0);
    double Z_pred     = (R0*Z0)(0,0);
    double gauss_Y    = Gauss(Y_pred, Y1(0,0), sigma_t);//Gauss2D(R0*Y0,Y1,Sigma); 
    double gauss_Z    = Gauss(Z_pred, Z1(0,0), sigma_t);//Gauss2D(R0*Z0,Z1,Sigma); 
    double gauss_WET  = Gauss(WEPL_mes, WEPL_est, sigma_WET);
    
    /*
    if(WEPL_mes>0 && WEPL_mes<2){
      cout<<" ------\nY0 "<<Y0(0,0)<<" -> "<<Y1(0,0)<<" Predicted :"<<(R0*Y0)(0,0)<<" PY0: "<<Y0(1,0)<<" PY1: "<<Y1(1,0)<<" Sigma_t: "<<sigma_t<<" Gauss_Y: "<<gauss_Y<<endl;
      cout<<"Z0 "<<Z0(0,0)<<" -> "<<Z1(0,0)<<" Predicted :"<<(R0*Z0)(0,0)<<" PZ0: "<<Z0(1,0)<<" PZ1: "<<Z1(1,0)<<" Sigma_t: "<<sigma_t<<" Gauss_Z: "<<gauss_Z<<endl;

      cout<<"P0 "<<p0.x()<<" "<<p0.y()<<" "<<p0.z()<<endl;
      cout<<"P1 "<<p1.x()<<" "<<p1.y()<<" "<<p1.z()<<endl;
      cout<<"Estop Est "<<Estop_est<<" Estop Mes: "<<Estop_mes<<" sigma_E: "<<sigma_E<<endl;
      cout<<"WEPL Est "<<WEPL_est<<" WEPL Meas "<<WEPL_mes<<" sigma_WET: "<<sigma_WET<<endl;
      cout<<"Prob Gauss: "<<gauss_Y<<" "<<gauss_Z<<" "<<gauss_WET<<endl;
      }
    */
    Point.WET_prob = gauss_WET;
    Point.Y_prob = gauss_Y;
    Point.Z_prob = gauss_Z;      
    bpt_WET->Fill();
    bpt_Y->Fill();
    bpt_Z->Fill();
  }

  t->Write("",TObject::kOverwrite);
  f->Close();
  
}

////////////////////////////////////////////
// Extract WET
////////////////////////////////////////////
double findWET(double Einit,double Estop){
  int it_Einit = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
  int it_Estop = lower_bound(Energy.begin(), Energy.end(), Estop)-Energy.begin();
  double WET = 0 ;
  for(int i=it_Estop;i<it_Einit;i++){
    WET += 0.01/dEdXBins[i]; // cm
  }
  return WET;
}
////////////////////////////////////////////
// Extract Energy
////////////////////////////////////////////
double findEnergy(double Einit,double WET){
  int it_Einit = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
  double temp_WET = 0;
  double Estop = 0;
  for(int i=it_Einit; i>0; i--){
    temp_WET += 0.01/dEdXBins[i]; // cm
    if( abs(WET-temp_WET) < 0.1 ) Estop = Energy[i];
  }
  return Estop;
}

////////////////////////////////////////////
// Compute Spline and Estimate Estop
////////////////////////////////////////////
double EstimateEstop(Particle *Point, TH3D* RSPMap, TMatrixD& vector_L, TMatrixD& vector_pv, TMatrixD& vector_E, TMatrixD& vector_X){
  TVector3 p0(Point->px0,Point->py0,Point->pz0);
  TVector3 p1(Point->px1,Point->py1,Point->pz1);
  TVector3 m0(Point->x0/10, Point->y0/10, Point->z0/10); // mm -> cm
  TVector3 m1(Point->x1/10, Point->y1/10, Point->z1/10); // mm -> cm

  TVector3 m, m_entry, m_exit, m_old; // position of the path at the entry and exit
  double t, t_entry = 0., t_exit = 1.; // fraction of the path at which point the spline enter the Hull

  // Negative rotation because we change the coordinates and not the phantom
  m0.RotateZ(-1*Point->angle*M_PI/180.);
  m1.RotateZ(-1*Point->angle*M_PI/180.);
  p0.RotateZ(-1*Point->angle*M_PI/180.);
  p1.RotateZ(-1*Point->angle*M_PI/180.);
  m_entry = m0;
  m_exit  = m1;
  
  //Problem from cubic root histogram and rotation here: rotating the object might leave the starting position out of the object. 
  //Hence the m0 must initially be forwared onto the reconstruction zylinder. 
  //The code assumes the center of the coordiante system lies in the center of the reconstructed volume
  //the reconstructed volume has a radius of radius_r = 12cm for the head, 9cm for Catphans
  double recon_radius = 15; //cm 15 for head, 9 for Catphan
 
  
  // Propagate from the entrance to the Hull
  for(int k =1; k<NStep-1; k++){
    t_entry = double(k)/NStep;
    m_entry = m0 + t_entry*p0;
    if (sqrt(pow(m_entry.x(),2)+pow(m_entry.y(),2))<recon_radius){ //Only now we are in the area covered by the reconstruction, nothing outside this radius
      int global = RSPMap->FindBin(m_entry.x(),m_entry.y(),m_entry.z());
      double RSP = RSPMap->GetBinContent(global);
      if(RSP>0.4) break;
    }
  }

  // Retro Propagate from the exit to the Hull
  for(int k =NStep-1; k>=1; k--){
    t_exit = double(k)/NStep;
    m_exit = m1 - t_exit*p1;
    if (sqrt(pow(m_exit.x(),2)+pow(m_exit.y(),2))<recon_radius){ //same as above
      int global = RSPMap->FindBin(m_exit.x(),m_exit.y(),m_exit.z());
      double RSP = RSPMap->GetBinContent(global);
      if(RSP>0.4) break;
    }
  }

  if(t_entry > t_exit) { // No Hull we are in air
    m_exit = m1;
    m_entry = m0;
    t_entry = 0.;
    t_exit = 1.;
  }

  double WER    = 25.693404; // 200 MeV -- cm
  double wepl   = Point->wepl/10; // mm -> cm
  double alpha1 = 1.01+0.43*pow(wepl/WER,2);
  double alpha2 = 0.99-0.46*pow(wepl/WER,2);
  double TotLength  = TVector3(m1-m0).Mag();
  double HullLength = TVector3(m_exit-m_entry).Mag(); 

  double RSP    = 0.0;
  double Einit  = 200;  
  p0.SetMag(alpha1*HullLength);
  p1.SetMag(alpha2*HullLength);
  TVector3 A,B,C,D;
  A       =    2*m_entry - 2*m_exit + p0+p1;
  B       =   -3*m_entry + 3*m_exit - 2*p0-p1;
  C       =    p0;
  D       =    m_entry;
  m_old   =    m0;
  for(int i=0;i<NStep;i++){
    
    t = double(i)/NStep;
    // before the Hull -- straight line
    if(t < t_entry) m = m0 + t*(m_entry-m0);
    
    // after the Hull -- straight line
    else if (t > t_exit) m = m_exit+ t*(m1-m_exit);

    // in the Hull -- cubic spline
    else m = D+t*(C+t*(B+t*A));

    // verifiy where we are in the prior reconstruction 
    if(m.x() < reco_Xmin || m.x() > reco_Xmax) RSP = 0.001; //RSP outside recon object not zero
    else if(m.y() < reco_Ymin || m.y() > reco_Ymax) RSP = 0.001;
    else if(m.z() < reco_Zmin || m.z() > reco_Zmax) RSP = 0.001;
    else{
      int binID  = RSPMap->FindBin(m.x(), m.y(), m.z()); 
      RSP = RSPMap->GetBinContent(binID);
    }      

    float L    = TVector3(m-m_old).Mag();   

    // fill the vectors to later calculate the sigmas
    if(RSP<0.4){ 
	vector_L(i,0) = 1./3.5E4; // cm
	RSP = 0.001;
    }
    else vector_L(i,0) = 1./(36.1); // cm
    vector_pv(i,0) = pow( (Einit+938.27)/((Einit + 1876.54)*Einit), 2);
    vector_E(i,0)  = Einit;    

    // to work with the rotation
    vector_X(i,0)  = t*TotLength;

    // keep a record of the exit energy
    int idE    = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
    Einit     -= (L)*RSP*dEdXBins[idE]*10;
    m_old      = m;
  }

  return Einit;
}
////////////////////////////////////////////
// Simpson rule of integration for the three sigma of y  along the x axis
////////////////////////////////////////////
void simps_mult(TMatrixD x, TMatrixD y, double t, double &sigma_t0, double &sigma_tt0, double &sigma_theta0){
  double h = 0;
  sigma_tt0 = 0; sigma_t0 = 0; sigma_theta0 = 0;
  for (int i=0; i < x.GetNoElements(); i++){

    if ( i == 0 || i == x.GetNoElements()-1 ){ // for the first and last elements
      h = (x(1,0)-x(0,0));
      sigma_theta0 += y(i,0)*h/3.;
      sigma_tt0    += (t-x(i,0))*y(i,0)*h/3.;
      sigma_t0     += (t-x(i,0))*(t-x(i,0))*y(i,0)*h/3.;
    }
    else
      {
	h = (x(i+1,0)-x(i,0));
	if (i%2==0)
	  {
	    sigma_theta0 += 2*y(i,0)*h/3.;
	    sigma_tt0    += 2*(t-x(i,0))*y(i,0)*h/3.;
	    sigma_t0     += 2*(t-x(i,0))*(t-x(i,0))*y(i,0)*h/3.;
	  }
	else
	  {
	    h = (x(i+1,0)-x(i,0));
	    sigma_theta0 += 4*y(i,0)*h/3.;
	    sigma_tt0    += 4*(t-x(i,0))*y(i,0)*h/3.;
	    sigma_t0     += 4*(t-x(i,0))*(t-x(i,0))*y(i,0)*h/3.;
	  }
      }

  }
}	
////////////////////////////////////////////
// element-wise multiplication of matrix
////////////////////////////////////////////
TMatrixD Mult(TMatrixD A, TMatrixD B){
  int nrows  = A.GetNrows();
  int ncols = A.GetNcols();
  TMatrixD result(nrows,ncols);
  for(int i=0;i<nrows ;i++){
    for(int j=0;j<ncols;j++){
      result(i,j) = A(i,j)*B(i,j);
    }
  }
  return result;
}
////////////////////////////////////////////
// element-wise multiplication
////////////////////////////////////////////
TMatrixD Mult(TMatrixD A, double B){
  int nrows  = A.GetNrows();
  int ncols = A.GetNcols();
  TMatrixD result(nrows,ncols);
  for(int i=0;i<nrows ;i++){
    for(int j=0;j<ncols;j++){
      result(i,j) = A(i,j)*B;
    }
  }
  return result;
}
////////////////////////////////////////////
// Simpson rule of integration
////////////////////////////////////////////
double simps(TMatrixD x, TMatrixD y){
  double sum=0.0;
  double h =0.0;
  for (int i=0; i < x.GetNoElements(); i+=1){

    if ( i == 0 || i == x.GetNoElements()-1 ) // for the first and last elements
      {
	h = x(1,0)- x(0,0);
	sum += y(i,0)*h/3.;
      }
    else
      {
	h = x(i+1,0)- x(i,0);
	if (i%2==0)
	  {
	    sum += 2*y(i,0)*h/3.;
	  }
	else
	  {
	    sum += 4*y(i,0)*h/3.; // the rest of data
	  }
      }
  }

  return sum;

}
////////////////////////////////////////////
// Energy straggling --> return the variance
////////////////////////////////////////////
double EnergyStraggling(TMatrixD tracks_E, double Eout){
  double strag = 0.;
  double K1     = 170./1000; // keV/cm -> MeV->cm
  double K2     = 0.082;//MeV2/cm
  double chi1, chi2, beta2;
  double dT;
  double mec2 = 0.511; // MeV
  double rho_e = 1.;
  double water_I = 78E-6;
  for(int idx=1;idx<NStep;idx++){
    dT     = abs(tracks_E(idx,0)- tracks_E(idx-1,0));
    beta2  = E2beta(tracks_E(idx,0));// beta2
    chi1   = (K1/beta2)*(log(2*mec2*beta2/((75E-6)*(1-beta2))) - beta2);
    chi2   = rho_e*K2*(1-0.5*beta2)/(1-beta2); // Tschalar
    strag += dT*chi2/(pow(chi1,3));
  }
  //Eout
  beta2    = E2beta(Eout);
  chi1     = (K1/beta2)*(log(2*mec2*beta2/((water_I)*(1-beta2))) - beta2);
  strag    = pow(chi1,2)*strag;
  return strag;
}
////////////////////////////////////////////
// Velocity divided by light velocity
////////////////////////////////////////////
double E2beta(double E){
  double mc2   = 938.27; // relativistic mass for protons
  double tau   = E/mc2;
  return (tau+2)*tau/pow(tau+1,2);

}
////////////////////////////////////////////
// Single-variate normal distribution
////////////////////////////////////////////
double Gauss(double x, double x0, double sigma){
  double num = TMath::Exp(-pow(x-x0,2)/(2*pow(sigma,2)));
  //double den = TMath::Sqrt(2*TMath::Pi()*pow(sigma,2));
  return  num;//den;
}
////////////////////////////////////////////
// Bi-variate normal distribution
////////////////////////////////////////////
double Gauss2D(TMatrixD Y0, TMatrixD Y1, TMatrixD Sigma){
  TMatrixD Sigma_I = Sigma.Invert();
  TMatrixD Diff  = Y1 - Y0;
  TMatrixD Diff_t(TMatrixD::kTransposed,Diff);  
  TMatrixD Part1 = Sigma_I*Diff;
  TMatrixD num   = (Diff_t*Part1);
  //double den     = (2*TMath::Pi())*TMath::Sqrt(Sigma.Determinant());
  return TMath::Exp(-0.5*num(0,0));///den;

}



