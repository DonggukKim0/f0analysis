#include <TFile.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TObject.h>
#include <TString.h>

void InvMass() {
    // TFile* Fin = new TFile("/Volumes/Seul/HP_results/242865.root","read"); 
    TFile* Fin = new TFile("/Users/seul/Downloads/AnalysisResults.root","read"); 

	THnSparse* hInvMassUS = (THnSparse*)Fin->Get("lf-f0980analysis/hInvMass_f0980_US");
	THnSparse* hInvMassLSpp = (THnSparse*)Fin->Get("lf-f0980analysis/hInvMass_f0980_LSpp");
	THnSparse* hInvMassLSmm = (THnSparse*)Fin->Get("lf-f0980analysis/hInvMass_f0980_LSmm");

    //  Multiplicity binning
    const int nmult = 8;
    int mult_min[nmult] = {0, 1, 5, 10, 20, 30, 50, 0};
	int mult_max[nmult] = {1, 5, 10, 20, 30, 50, 100, 100};
    //  pT binning
    const int npt = 16;
    double pt_min[npt] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
	double pt_max[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

    TH1D* hProjInvMassUS[nmult][npt];
	TH1D* hProjInvMassLSpp[nmult][npt];
	TH1D* hProjInvMassLSmm[nmult][npt];
	TH1D* hProjInvMassLS[nmult][npt];
	TH1D* hProjInvMassSub[nmult][npt];

    for (int i = 7; i < nmult; i++) {
        hInvMassUS->GetAxis(2)->SetRangeUser(mult_min[i], mult_max[i]-0.001);
		hInvMassLSpp->GetAxis(2)->SetRangeUser(mult_min[i], mult_max[i]-0.001);
		hInvMassLSmm->GetAxis(2)->SetRangeUser(mult_min[i], mult_max[i]-0.001);

        TCanvas *c = new TCanvas(Form("c_%d", i), Form("Canvas %d", i), 1500, 700);
		c->Divide(6,3,0.001,0.001);
		c->SetGrid();

        for (int j = 0; j < npt; j++) {
            hInvMassUS->GetAxis(1)->SetRangeUser(pt_min[j], pt_max[j]-0.001);
			hInvMassLSpp->GetAxis(1)->SetRangeUser(pt_min[j], pt_max[j]-0.001);
			hInvMassLSmm->GetAxis(1)->SetRangeUser(pt_min[j], pt_max[j]-0.001);

            //  Get the mass distribution - 0 axis, contain erorr "e"
            hProjInvMassUS[i][j] = hInvMassUS->Projection(0,"e");
			hProjInvMassUS[i][j]->SetName(Form("hInvMassUS_%d_%d", i, j));
			hProjInvMassLSpp[i][j] = hInvMassLSpp->Projection(0,"e");
			hProjInvMassLSpp[i][j]->SetName(Form("hProjInvMassLSpp_%d_%d",i, j));
			hProjInvMassLSmm[i][j] = hInvMassLSmm->Projection(0,"e");
			hProjInvMassLSmm[i][j]->SetName(Form("hProjInvMassLSmm_%d_%d",i, j));

			hProjInvMassLS[i][j] = (TH1D*)hProjInvMassUS[i][j]->Clone();
			hProjInvMassLS[i][j]->SetName(Form("hProjInvMassLS_%d_%d",i, j));
			hProjInvMassLS[i][j]->Reset();
            //  LS fill - 2*sqrt(Npp+Nmm)
			for (int p = 0; p < hProjInvMassUS[i][j]->GetNbinsX(); p++) {
                double LSpp = hProjInvMassLSpp[i][j]->GetBinContent(p+1);
                double LSmm = hProjInvMassLSmm[i][j]->GetBinContent(p+1);
                double LSpp_err = hProjInvMassLSpp[i][j]->GetBinError(p+1)/LSpp;
                double LSmm_err = hProjInvMassLSmm[i][j]->GetBinError(p+1)/LSmm;

				hProjInvMassLS[i][j]->SetBinContent(p+1, 2.0*sqrt(LSpp*LSmm));
				hProjInvMassLS[i][j]->SetBinError(p+1, sqrt(LSpp*LSmm)*(pow(LSpp_err,2)+pow(LSmm_err,2)));
			}
            //  US, LS subtraction
			hProjInvMassSub[i][j] = (TH1D*)hProjInvMassUS[i][j]->Clone();
			hProjInvMassSub[i][j]->GetXaxis()->SetRangeUser(0.5, 1.75);
			hProjInvMassSub[i][j]->SetName(Form("hProjInvMassSub_%d_%d", i, j));
            hProjInvMassSub[i][j]->Add( hProjInvMassLS[i][j], -1.0 );
   
            //  Draw
            const int cIdx = j+1;
            c->cd(cIdx);
            gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);

			hProjInvMassSub[i][j]->Draw("hist");
			// hProjInvMassSub[i][j]->Rebin(2); // Binning 2000 -> 1000
			hProjInvMassSub[i][j]->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
			// hProjInvMassSub[i][j]->GetXaxis()->SetRangeUser(0.8, 1.75);
			hProjInvMassSub[i][j]->GetXaxis()->SetRangeUser(0.7, 1.75);

			TLegend* leg = new TLegend(0.4,0.6,0.9,0.85); // x1, y1, x2, y2
			leg->SetTextSize(0.045);
			leg->SetLineWidth(0.0);
			leg->SetFillStyle(0);
			leg->AddEntry((TObject*)0, "LHC22o apass7", "");
			leg->AddEntry((TObject*)0, Form("FT0M %d#font[122]{-}%d %%", mult_min[i], mult_max[i]), "");
			leg->AddEntry((TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.5", "");
			// leg->AddEntry( (TObject*)0, Form("#it{p}_{T,lead} > 5 GeV/#it{c}, %s",trnsName[r]), "");
			leg->AddEntry((TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",pt_min[j], pt_max[j]), "");
			leg->Draw();
        }
        TString fileName = Form("plot/Invmass_mult_%d_%d.pdf", mult_min[i], mult_max[i]);
		c->SaveAs(fileName);
		delete c;
    }

    // Create root output file
    TFile* Fout = new TFile("/Users/seul/seul_workspace/f0ana/results/InvMassOut.root","recreate");
    for (int i = 7; i < nmult; i++) {
        for (int j = 0; j < npt; j++) {
            hProjInvMassUS[i][j]->Write();
            hProjInvMassLS[i][j]->Write();
            hProjInvMassSub[i][j]->Write();
        }
    }
}
