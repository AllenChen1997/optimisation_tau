/////////////////////////////////////////////////
/*
   this version only do the n2b1 first
*/
////////////////////////////////////////////////

#include "tau21_study.C"

using namespace std;

string filesdir = "/afs/cern.ch/work/k/kuchen/public/samples/tauver/";
string sig_root= filesdir+"crab_EXO-ggToXdXdHToBB_sinp_0p35_tanb_1p0_mXd_10_MH3_1000_MH4_150_MH2_1000_MHC_1000_CP3Tune_13TeV_0000_1.root";
string tt_semi= filesdir+"crab_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.root";
string tt_LL = filesdir+"crab_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.root";
string tt_had = filesdir+"crab_TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8.root";

string s1 = "h_pt_rho_26";
bool isddt = false;
double Mintau21;

// the top sample xs:
#define semi 687.1*0.438
#define LL 687.1*0.105
#define hadron 687.1*0.457
#define Lumi 41000.0 // 1/pb

void load_to_hist(string s, TH1D* h, TH2D* h_cut, double xsbkg, bool isQCD){
	TFile* myfile = new TFile(s.c_str(),"READ");
	TTreeReader myRead("monoHbb_SR_boosted",myfile);  
	TTreeReaderValue< Double_t > tau21(myRead,"fjetTau21");
	TTreeReaderValue< Double_t > rho(myRead,"FJetrho");
	TTreeReaderValue< Double_t > pt(myRead,"FJetPt");
	TTreeReaderValue< Double_t > dphi(myRead,"min_dPhi");
	TTreeReaderValue< Double_t > ddb(myRead,"FJetCSV");
	TTreeReaderValue< Double_t > nj(myRead,"nJets");
	TTreeReaderValue< Double_t > mass(myRead,"FJetMass");
	double width = (double) (Maxrho -Minrho) / (double) Nrho;
	while (myRead.Next()){  // loop in one root file
		/*if (*dphi < 0.4 ) continue;
		if (*mass < 100 || *mass > 150) continue;
		if (*ddb < 0.86 ) continue;
		if (*nj > 2) continue;*/
		if (isddt){
			if (*rho < Minrho ||*rho > Maxrho) continue;
			int x = ceil((double)(*rho - Minrho) / d );
			int y;
			if (*pt >= pt_r1 && *pt < pt_r2 ) y = 1;
			else if (*pt >= pt_r2 && *pt < pt_r3) y = 2;
			else if (*pt >= pt_r3 && *pt < pt_r4) y = 3;
			else if (*pt >= pt_r4) y = 4;
			else continue;
			double n_cut = h_cut->GetBinContent(x,y);
			if (n_cut == 0) continue;
			h->Fill(*tau21-n_cut);
		} else {
			h->Fill(*tau21);
		}
	}
	
	if (xsbkg != 0){	
		TH1F* h_event = (TH1F*) myfile->Get("h_total_mcweight");
		double totalevent = h_event->Integral();
		if (! isQCD ) h->Scale(1.0/totalevent);
	}
}
void plot_3(){
	if (isddt) Mintau21 = -0.5;
	else Mintau21 = 0;
	TH1D* h_sig = new TH1D("signal","signal",Ntau,Mintau21,Maxtau);
	TH1D* h_top[3];
	h_top[0] = new TH1D("top_1","top_1",Ntau,Mintau21,Maxtau);
	h_top[1] = (TH1D*) h_top[0]->Clone("top_2");
	h_top[2] = (TH1D*) h_top[0]->Clone("top_3");
	TH1D* h_QCD = new TH1D("QCD","QCD",Ntau,Mintau21,Maxtau);
	
	TFile* myfile = new TFile("TH3_output.root","READ");
	TH2D* h_pt_rho = (TH2D*) myfile->Get(s1.c_str());
	
	// QCD //
	ifstream infile("QCD_xs_list.txt");
	string line,name;
	double xs;
	stringstream ss;
	while(getline(infile,line)){
		ss << line;
		ss >> name >> xs;
		ss.clear();
		load_to_hist(name,h_QCD,h_pt_rho,xs,true);
	}
	// signal //
	load_to_hist(sig_root,h_sig,h_pt_rho,0,false); // xs = 0, used to show it is signal in code
	// top //
	double eff=0;
	load_to_hist(tt_semi,h_top[0],h_pt_rho,semi,false);
	load_to_hist(tt_LL,h_top[1],h_pt_rho,LL,false);
	load_to_hist(tt_had,h_top[2],h_pt_rho,hadron,false);
	// output plot //
	h_top[0]->Scale(semi);
	h_top[1]->Scale(LL);
	h_top[2]->Scale(hadron);

	gStyle->SetOptStat("");	
	h_top[0]->Add(h_top[1]);
	h_top[0]->Add(h_top[2]);
	
	auto c1 = new TCanvas("c1","c1");
	h_sig->Scale(1.0/h_sig->Integral());
	if (isddt) h_sig->SetXTitle("tau21DDT");
	else h_sig->SetXTitle("tau21");
	h_sig->SetTitle("");
	h_sig->SetLineColor(kRed);
	h_sig->GetXaxis()->SetTitleSize(0.04);
	h_sig->GetXaxis()->SetLabelSize(0.05);
	h_sig->GetYaxis()->SetLabelSize(0.05);
	h_sig->SetTitleOffset(1.0,"X");
	
	h_sig->Draw("HIST E");
	h_QCD->Scale(1.0/h_QCD->Integral());
	h_QCD->Draw("SAME HIST E");
	h_top[0]->Scale(1.0/h_top[0]->Integral());
	h_top[0]->SetLineColor(kBlack);
	h_top[0]->Draw("SAME HIST E");
	
	TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
	legend->AddEntry(h_sig,"signal","l");
	legend->AddEntry(h_QCD,"QCD bkg","l");
	legend->AddEntry(h_top[0],"top bkg","l");
	legend->Draw();
	
	if (isddt) c1->SaveAs("tau21DDT_new2.png");
	else c1->SaveAs("tau21_new2.png");
}
