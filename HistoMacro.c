void HistoMacro(){
   TFile *f = new TFile("output.root");
   TCanvas *c = new TCanvas("c","c", 800, 700);

   gStyle->SetOptStat(0);
   //gStyle->SetOptTitle(0);
   gStyle->SetFrameLineWidth(2);
   gStyle->SetLabelOffset(1.2);
   gStyle->SetLabelFont(72);

//Histograms
   TH1F *h1Z = (TH1F*)f->Get("Zmass_MllCut_WZQuad_all;1");
   TH1F *h2Z = (TH1F*)f->Get("Zmass_MllCut_SMWZqcd_all;1");
   TH1F *h3Z = (TH1F*)f->Get("Zmass_MllCut_SMWZew_all;1");

   TH1F *h1W = (TH1F*)f->Get("Wmass_MllCut_WZQuad_all;1");
   TH1F *h2W = (TH1F*)f->Get("Wmass_MllCut_SMWZqcd_all;1");
   TH1F *h3W = (TH1F*)f->Get("Wmass_MllCut_SMWZew_all;1");

   TH1F *h1WZ = (TH1F*)f->Get("WZmass_MllCut_WZQuad_all;1");
   TH1F *h2WZ = (TH1F*)f->Get("WZmass_MllCut_SMWZqcd_all;1");
   TH1F *h3WZ = (TH1F*)f->Get("WZmass_MllCut_SMWZew_all;1");

   //h1Z->Scale(1./h1Z->Integral()); // normalization
   //h2Z->Scale(1./h2Z->Integral()); // normalization
   //h3Z->Scale(1./h3Z->Integral()); // normalization  


   //histogram editing options :
   h1Z->SetFillStyle(1001);
   h1Z->SetFillColorAlpha(kRed+1,0.33);
   h1Z->SetTitle("Zmass");
   h1Z->GetXaxis()->SetTitle("GeV");
   h1Z->GetYaxis()->SetTitle("");
   h1Z->GetYaxis()->SetRangeUser(0,5 * pow(10,-2));
   h1Z->GetXaxis()->SetRangeUser(0,0);
   h1Z->SetLineWidth(2);
   h1Z->SetLineColor(kRed+1);

   h2Z->SetFillStyle(1001);
   h2Z->SetFillColorAlpha(kBlue+1,0.33);
   h2Z->SetLineWidth(2);
   h2Z->SetLineColor(kBlue+1); 
   
   h3Z->SetFillStyle(1001);
   h3Z->SetFillColorAlpha(kGreen+1,0.33);
   h3Z->SetLineWidth(2);
   h3Z->SetLineColor(kGreen+1);  

   h1Z->Draw("Hist");
   h2Z->Draw("Hist" "same");
   h3Z->Draw("Hist" "same");

   TLatex *t = new TLatex();
   t->SetTextFont(42);
   t->SetTextSize(0.028);
   //t->DrawLatex(7.840,0.27," text_if_needed ");
   t->Draw("same");
   
   auto legend = new TLegend(0.666,0.704,0.840,0.838);
   legend->AddEntry(h1Z,"WZ Quad","f");
   legend->AddEntry(h2Z,"WZ SM qcd","f");
   legend->AddEntry(h3Z,"WZ SM ew","f");
   gStyle->SetLegendFont(42);
   gStyle->SetLegendTextSize(0.03);
   gStyle->SetLegendBorderSize(0);
   legend->SetFillStyle(0);
   legend->SetBorderSize(0);

   legend->Draw("");
   c->Print("Histogram.pdf["); // Opens pdf
   c->Print("Histogram.pdf"); //prints the first page
   c->Clear();

   h1W->SetFillStyle(1001);
   h1W->SetFillColorAlpha(kRed+1,0.33);
   h1W->GetXaxis()->SetTitle("GeV");
   h1W->GetYaxis()->SetTitle("");
   h1W->GetYaxis()->SetRangeUser(0,0.012);
   h1W->GetXaxis()->SetRangeUser(0,0);
   h1W->SetLineWidth(2);
   h1W->SetLineColor(kRed+1); 

   h2W->SetFillStyle(1001);
   h2W->SetFillColorAlpha(kBlue+1,0.33);
   h2W->SetLineWidth(2);
   h2W->SetLineColor(kBlue+1);

   h3W->SetFillStyle(1001);
   h3W->SetFillColorAlpha(kGreen+1,0.33);
   h3W->SetLineWidth(2);
   h3W->SetLineColor(kGreen+1);  

   h1W->Draw("Hist");
   h2W->Draw("Hist" "Same");
   h3W->Draw("Hist" "Same");

   t->Draw("same");
   legend->Draw("");
   c->Print("Histogram.pdf");

   h1WZ->SetFillStyle(1001);
   h1WZ->SetFillColorAlpha(kRed+1,0.33);
   h1WZ->GetXaxis()->SetTitle("GeV");
   h1WZ->GetYaxis()->SetTitle("");
   h1WZ->GetYaxis()->SetRangeUser(0,0.012);
   h1WZ->GetXaxis()->SetRangeUser(0,0);
   h1WZ->SetLineWidth(2);
   h1WZ->SetLineColor(kRed+1); 

   h2WZ->SetFillStyle(1001);
   h2WZ->SetFillColorAlpha(kBlue+1,0.33);
   h2WZ->SetLineWidth(2);
   h2WZ->SetLineColor(kBlue+1);

   h3WZ->SetFillStyle(1001);
   h3WZ->SetFillColorAlpha(kGreen+1,0.33);
   h3WZ->SetLineWidth(2);
   h3WZ->SetLineColor(kGreen+1);  

   h1WZ->Draw("Hist");
   h2WZ->Draw("Hist" "Same");
   h3WZ->Draw("Hist" "Same");

   t->Draw("same");
   legend->Draw("");
   c->Print("Histogram.pdf");
 
   c->Print("Histogram.pdf]"); //Closes pdf




}