void MergedHists()
{
const int num = 4;
vector<Color_t> colours = {kBlue, kRed, kGreen, kOrange};
const char* grtitles[num] = {"Btagging+MVA", "NoBtagging+MVA", "BtaggingNoMVA", "NoBtaggingNoMVA"};
const char* dirnames[num] = {"LeadingCut/","NoLeadingCut/", "LeadingNoCut/", "NoLeadingNoCut/"};
const char* sigfiles[num] = {"TMVApp_sigMore25.root","TMVApp_sig2550.root", "TMVApp_sig5075.root", "TMVApp_sigMore75.root" };
const char* bkgfiles[num] = {"TMVApp_bkgMore25.root", "TMVApp_bkg2550.root",  "TMVApp_bkg5075.root", "TMVApp_bkgMore75.root"};
const char* ROC[num] = {"ROCMore25", "ROC2550",  "ROC5075", "ROCMore75"};
const char* Acc[num] = {"AccMore25", "Acc2550",  "Acc5075", "AccMore75"};
auto mg = new TMultiGraph();
TGraph* gr[num];
const int n=100;
Double_t x[n], y[n];

    TCanvas* c1 = new TCanvas();
    gStyle->SetOptStat(0);
for (int i=0; i<num; i++)
{
char title1[100], title2[100];
strcpy(title1, dirnames[i]);
strcat(title1, sigfiles[3]);

strcpy(title2, dirnames[i]);
strcat(title2, bkgfiles[3]);
cout<<title1<<" "<<title2<<endl;


    TFile *fs = new TFile(title1, "READ");
    TFile *fb = new TFile(title2, "READ");
    TH1D * hb = new TH1D("hb","background", 20, 0, 1.5);
    TH1D * hs = new TH1D("hs","signal", 20, 0, 1.5);

    hs = (TH1D*)fs->Get("MVA_BDT");
    hb = (TH1D*)fb->Get("MVA_BDT");

    TH1* h1 = hs->GetCumulative();
    TH1* h2 = hb->GetCumulative();
    h1->SetMarkerStyle(kFullCircle);
    h1->SetMarkerSize(1);
    h2->SetMarkerSize(1);
    h2->SetMarkerStyle(kFullSquare);
    h1->SetMarkerColor(kBlue);
    h2->SetMarkerColor(kGreen);
    h1->Scale(1/hs->Integral());
    h2->Scale(1/hb->Integral());
/*    h1->Draw("P");
    h2->Draw("SAME P");
    c1->Print("Acc.png");
    c1->Print("Acc.C");*/
    for (int j=1; j<n+1; j++)
    {
	    x[j] = 1 - h1->GetBinContent(j);
	    y[j] = h2->GetBinContent(j);
	    cout<<x[j]<<" "<<y[j]<<" "<<h1->GetBinCenter(i)<<" "<<h2->GetBinCenter(j)<<endl;
    }
	gr[i] = new TGraph(n, x, y);
	gr[i]->SetLineColor(colours[i]);
	gr[i]->SetMarkerColor(colours[i]);
	gr[i]->SetMarkerSize(1);
	gr[i]->Draw("AC*");

        mg->Add(gr[i]);

	delete hs;
	delete hb;
	delete h1;
	delete h2;
	delete fs;
	delete fb;
}
//mg->SetName("More25");
mg->GetXaxis()->SetTitle("Signal");
mg->GetYaxis()->SetTitle("Rejected background");
mg->Draw("AP");
c1->Update();
//c1->BuildLegend();
   c1->Print("ROCMore75.png");
   c1->Print("ROC.C");

}
