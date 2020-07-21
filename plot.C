///////////////////////////////
/*
	Date 2020/ 7/ 21
	owner Kung-Hsiang Chen
	
	read signal and QCD samples and plot tau21 figures
*/
///////////////////////////////
using namespace std;
#define Ntau 10
#define Mintau 0.0
#define Maxtau 1.0
#define NN 8 // the number of histos.

template<typename T>
void setDrawOpt(T& h,string title, string xTitle, string yTitle){
	h->SetTitle(title.c_str());
	h->SetTitleSize(0.07);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->SetXTitle(xTitle.c_str());
	h->SetYTitle(yTitle.c_str());
}

void getPlot(string inputtxt, string outputname){
	TH1D* h_out = new TH1D("hout","",Ntau,Mintau,Maxtau); setDrawOpt(h_out,"","#tau","");
	TH1D* h_tmp[NN];
	for (int i=0;i<NN;i++) h_tmp[i] = (TH1D*) h_out->Clone("tmp"); 
	ifstream infile(inputtxt); // used to input file. in each line: xxx.root name_for_plot
	string line,s1,s2;
	stringstream ss;
	auto c1 = new TCanvas("c1","c1");
	TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
	int i=0;
	while(getline(infile,line)){
		cout << line << endl;
		ss << line;
		ss >> s1 >> s2 ;// s1=xxx.root, s2=name_for_plot
		ss.clear();
		TFile* myfile = new TFile(s1.c_str(),"READ");
		TTreeReader myRead("monoHbb_SR_boosted",myfile);  
		TTreeReaderValue< Double_t > tau(myRead,"fjetTau21");
		TTreeReaderValue< Double_t > rho(myRead,"FJetrho");
		TTreeReaderValue< Double_t > pt(myRead,"Jet1Pt");

		int N = myRead.GetEntries(); //get the entires info.
		if (N == 0) continue;
		while (myRead.Next()){
			h_tmp[i]->Fill(*tau);
		}
		// plot out the plots
		h_tmp[i]->Scale(1.0/(float)h_tmp[i]->Integral());
		h_tmp[i]->SetLineColor(i+1);
		h_tmp[i]->Draw("SAME");
		legend->AddEntry(h_tmp[i],(TString)s2,"l");
		i++;
		if (i==NN) {
			cout << "over the maximum histo. limit, please check NN in code" << endl;
			break;
		}
	}
	legend->Draw();
	c1->SaveAs((TString)outputname+".png");
	c1->Close();
	
}

void plot(){
	TH1D* h_signal = new TH1D("signal","",Ntau,Mintau,Maxtau);
	TH1D* h_QCD = (TH1D*) h_signal->Clone("QCD");
	gStyle->SetOptStat("");
	getPlot("./signal_list_for_plot.txt","signal_tau");
	getPlot("./QCD_list_for_plot.txt","QCD_tau");
}
