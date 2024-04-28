void HistoMacro(){
   TFile *f = new TFile("output.root");
   TCanvas *c = new TCanvas("c","c", 800, 700);

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetFrameLineWidth(2);
   gStyle->SetLabelOffset(1.2);
   gStyle->SetLabelFont(72);

   TH1F *h1 = (TH1F*)f->Get("Zmass_MllCut_WZQuad_all;1");
   TH1F *h2 = (TH1F*)f->Get("Zmass_MllCut_SMWZqcd_all;1");
   TH1F *h3 = (TH1F*)f->Get("Zmass_MllCut_SMWZew_all;1");

   h1->Scale(1./h1->Integral()); // normalization
   h2->Scale(1./h2->Integral()); // normalization
   h3->Scale(1./h3->Integral()); // normalization  


   //histogram editing options :
   h1->SetFillStyle(1001);
   h1->SetFillColorAlpha(kRed+1,0.33);
   h1->GetXaxis()->SetTitle("MeV");
   h1->GetYaxis()->SetTitle("");
   h1->GetYaxis()->SetRangeUser(0,1);
   h1->GetXaxis()->SetRangeUser(0,0); 

   h2->SetFillStyle(1001);
   h2->SetFillColorAlpha(kBlue+1,0.33);
   
   h3->SetFillStyle(1001);
   h3->SetFillColorAlpha(kGreen+1,0.33);

   h1->SetLineWidth(2);
   h1->SetLineColor(kRed+1);

   h2->SetLineWidth(2);
   h2->SetLineColor(kBlue+1); 

   h3->SetLineWidth(2);
   h3->SetLineColor(kGreen+1);   
   h1->Draw("Hist");
   h2->Draw("Hist" "same");
   h3->Draw("Hist" "same");



   TLatex *t = new TLatex();
   t->SetTextFont(42);
   t->SetTextSize(0.028);
   //t->DrawLatex(7.840,0.27," text_if_needed ");
   t->Draw("same");
   
   auto legend = new TLegend(0.666,0.704,0.840,0.838);
   legend->AddEntry(h1,"WZ Quad","f");
   legend->AddEntry(h2,"WZ SM qcd","f");
   legend->AddEntry(h3,"WZ SM ew","f");
   gStyle->SetLegendFont(42);
   gStyle->SetLegendTextSize(0.03);
   gStyle->SetLegendBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetBorderSize(0);

   legend->Draw("");
   c->Print("Zmass_Histogram.pdf");
}