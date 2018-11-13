#include <iostream>
#include "THStack.h"
#include <TMath.h>
#include "TLegend.h"
#include "xcdf/utility/NumericalExpression.h"
#include "TGraph.h"
#include "xcdf/XCDFField.h"
#include "TTree.h"
#include "xcdf/XCDF.h"
#include "TFile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TRint.h"
#include "vector"
#include "utility"
#include "TH2F.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TPaveLabel.h"
using namespace std;

const double hHAWC = 640.; // g.cm-2
const vector<char> tID ={'A','B','C','D'};
const float r2d = 360/(2*3.141592654);
const float d2r = 1/r2d;
const vector<int> col ={kBlack,kGreen,kBlue,kRed,kOrange,kGray+1};
//const Float_t fbins[] = {.485, .618, .74};
const Float_t fbins[] = {.067, .105, .162, .247, .356, .485, .618, .74, .84, 1};
const Float_t ebins[] = {.25, .5, .75, 1.0, 1.25, 1.5, 1.75};
//const Float_t ebin[] = {.25, .5, .75, 1.0, 1.25, 1.5, 1.75};
const Int_t  fbinnum = sizeof(fbins)/sizeof(Float_t) - 1; 
const Int_t  ebinnum = sizeof(ebins)/sizeof(Float_t) - 1; 

void getx(string name, TH1F** h1, int k, int end) {
    XCDFFile recFile(name.c_str(),"r");
    
	for(float i=0; i<=7; i+=1){
		float n=i/4;
		string b = to_string(n+.25);
		string nam = "Ebin " + b.substr(0,4);
		h1[(int) i] = new TH1F(nam.c_str(),nam.c_str(),25,0,end);
					}
						
	// First way to read data, a la SteBranchAddress
	XCDFUnsignedIntegerField nHit    	= recFile.GetUnsignedIntegerField("rec.nHitSP20");
	XCDFUnsignedIntegerField nch    	= recFile.GetUnsignedIntegerField("rec.nChAvail");
	XCDFUnsignedIntegerField fidu    	= recFile.GetUnsignedIntegerField("rec.coreFiduScale");
	XCDFFloatingPointField   zenAngle 	= recFile.GetFloatingPointField("rec.zenithAngle");
	XCDFFloatingPointField   iw 		= recFile.GetFloatingPointField("sweets.IWgt");
	//XCDFFloatingPointField   charge   = recFile.GetFloatingPointField("event.hit.charge");
	//XCDFFloatingPointField   effcharge= recFile.GetFloatingPointField("event.hit.charge");
	XCDFFloatingPointField 	 E          = recFile.GetFloatingPointField("mc.logEnergy");
	XCDFFloatingPointField   cxpe 		= recFile.GetFloatingPointField("rec.CxPE40");
	
		while (recFile.Read()) {
			double fhit = (double) (*nHit)/(*nch);
			double Et = *E-3;
			double v = *zenAngle;
			int fid = *fidu;
			
				if((fhit>fbins[k+1]) and (fhit<fbins[k+2])){
						for (int i=0;i<=6;i++){
							if((Et>ebins[i]) and (Et<ebins[i+1])){
								h1[i]->Fill(v,*iw);
																	}
											}
														}
				
			
			

								}//read
	
	}		

void norm(TH1F* h){
	double t= h->Integral();
	h->Scale(1.0/t);
	}
	
int main(){


	string name = "Zenith angle distribution between fhit between 0.485 and 0.618 over various Log(E) bins";
	TApplication app("app",0,0);
	TCanvas *c1 = new TCanvas();
	/* float xh = 2.5;
	float xl = 0.5;
	int binN = 40;*/
	float h = .225; 
	TH1F *hist[ebinnum+1];
	TGraph *grph[ebinnum+1];
	string xlab="Zenith Angle";
	string ylab="Frequency";
	TMultiGraph *mg = new TMultiGraph();
	auto leg = new TLegend(0.65, 0.65, .9, .9);
	

	int n = ebinnum;
	Double_t x[ebinnum+1][fbinnum+1], y[ebinnum+1][fbinnum+1];
	THStack *hs = new THStack("hs","");
	int k = 5;
	getx("g-out-file.xcd",hist,k,1);
	int j=k-4;
	for(int i = j; i<=ebinnum-1; i++){
		norm(hist[i]);
		hist[i]->GetYaxis()->SetTitle(ylab.c_str());
		hist[i]->GetXaxis()->SetTitle(xlab.c_str());
		hist[i]->SetLineWidth(2);
		hist[i]->SetLineColor(col[i]);
		hist[i]->GetYaxis()->SetRangeUser(0,h);
		c1->cd(i);
		
		cout<<"hist["<<i<<"]: ";
		for(int k = 0;k<=hist[i]->GetSize();k++){
			cout<<hist[i]->GetBinContent(k)<<", ";
			}
			cout<<endl;
		if(i == j){
			hist[i]->Draw();
			}
		else{
		hist[i]->Draw("sames");
		//hs->Add(hist[i]);
		} 
	}
	
	TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,name.c_str(),"brndc");
	title->Draw(); 

	//mg->Draw("ACP");
	//hs->Draw("nostack");

	c1->BuildLegend();
	c1->Paint();
	c1->Modified();
	c1->cd();
	c1->Update();
	app.Run();
	}

