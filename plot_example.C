#include "plot_helpers.h"


void plot_example() {
    // NOTE: These two lines are only needed if the macro is executed
    // outside of the build directory
    gSystem->AddDynamicPath("build/");
    gSystem->Load("libG__ClassesDict.so");

    const auto c = 299792458.; // m/s
    const auto pi = 3.14159265359;

    auto f = TFile::Open("out.root");

    auto tree = (TTree*)f->Get("nominal");

    std::vector<int>* prod_vtx  = 0;
    std::vector<int>* decay_vtx = 0;

    std::vector<std::vector<int>>* children = 0;
    std::vector<std::vector<int>>* parents  = 0;

    std::vector<double>* m        = 0;
    std::vector<double>* pt       = 0;
    std::vector<double>* e        = 0;
    std::vector<double>* eta      = 0;
    std::vector<double>* phi      = 0;
    std::vector<int>* pdg_id      = 0;
    std::vector<int>* vtx_barcode = 0;
    std::vector<double>* vtx_x    = 0;
    std::vector<double>* vtx_y    = 0;
    std::vector<double>* vtx_z    = 0;

    tree->SetBranchAddress("m", &m);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("e", &e);
    tree->SetBranchAddress("eta", &eta);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("prod_vtx", &prod_vtx);
    tree->SetBranchAddress("decay_vtx", &decay_vtx);
    tree->SetBranchAddress("pdg_id", &pdg_id);
    tree->SetBranchAddress("vtx_x", &vtx_x);
    tree->SetBranchAddress("vtx_y", &vtx_y);
    tree->SetBranchAddress("vtx_z", &vtx_z);
    tree->SetBranchAddress("vtx_barcode", &vtx_barcode);

    tree->SetBranchAddress("children", &children);
    tree->SetBranchAddress("parents", &parents);

    auto h_m        = new TH1D("h_mass", "", 100, 0, 500);
    auto h_e        = new TH1D("h_e", "", 100, 0, 500);
    auto h_m_mumu   = new TH1D("h_mumu_mass", "", 100, 50, 200);
    auto h_pt       = new TH1D("h_pt", "", 100, 0, 500);
    auto h_pt0      = new TH1D("h_pt0", "", 100, 0, 400);
    auto h_pt1      = new TH1D("h_pt1", "", 100, 0, 400);
    auto h_pt2      = new TH1D("h_pt2", "", 100, 0, 400);
    auto h_dl       = new TH1D("h_dl", "", 200, 0, 4000);
    auto h_dl0      = new TH1D("h_dl0", "", 200, 0, 1000);
    auto h_vtx_x    = new TH1D("h_vtx_x", "", 100, 0, 10000);
    auto h_gamma    = new TH1D("h_gamma", "", 100, 0, 50);
    auto h_children = new TH1D("h_children", "", 220, 0, 220);
    auto h_dr12     = new TH1D("h_dr12", "", 100, 0, 4);
    auto h_dr23     = new TH1D("h_dr23", "", 100, 0, 4);
    auto h_dr13     = new TH1D("h_dr13", "", 100, 0, 4);
    auto h_eta      = new TH1D("h_eta", "", 100, -4, 4);
    auto h_phi      = new TH1D("h_phi", "", 100, -pi, pi);

    auto h_eta_phi = new TH2D("h_eta_phi", "", 100, -4, 4, 100, -pi, pi);

    auto nevents = tree->GetEntries();
    for (auto i = 0; i < nevents; ++i) {
        tree->GetEntry(i);

        if (i % 500 == 0) {
            std::cout << "Event " << i << '\n';
        }

        auto nparts = pdg_id->size();
        for (auto ipart = 0; ipart < nparts; ++ipart) {
            auto id = (*pdg_id)[ipart];

            if (abs(id) == 23) {
                for (auto child : (*children)[ipart]) {
                    auto tmp = (*pdg_id)[child];
                    h_children->Fill(abs(tmp));
                }
            }

            if (abs(id) == 23 && is_first_parent(ipart, *parents, *pdg_id)) {
                auto jpart = find_last_child(ipart, *children, *pdg_id);

                auto v = to_v4(jpart, *pt, *eta, *phi, *m);

                h_pt->Fill(v.Pt() / 1000.);
                h_m->Fill(v.M() / 1000.);
                h_eta->Fill(v.Eta());
                h_phi->Fill(v.Phi());
                h_eta_phi->Fill(v.Eta(), v.Phi());
                h_m->Fill(v.M() / 1000.);
                h_e->Fill(v.E() / 1000.);

                double gamma = v.Gamma();
                h_gamma->Fill(gamma);

                auto vtx0 = (*prod_vtx)[ipart];
                auto vtx1 = (*decay_vtx)[jpart];

                double x0 = (*vtx_x)[vtx0];
                double y0 = (*vtx_y)[vtx0];
                double z0 = (*vtx_z)[vtx0];

                double x1 = (*vtx_x)[vtx1];
                double y1 = (*vtx_y)[vtx1];
                double z1 = (*vtx_z)[vtx1];

                double dl =
                        std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2) +
                                  std::pow(z1 - z0, 2));
                double dl0 = dl / gamma;

                h_dl->Fill(dl);
                h_dl0->Fill(dl0);

                if (children[jpart].size() > 1) {
                    std::vector<int> child = {
                            (*children)[jpart][0], (*children)[jpart][1],
                    };

                    auto v0 = to_v4(child[0], *pt, *eta, *phi, *m);
                    auto v1 = to_v4(child[1], *pt, *eta, *phi, *m);

                    h_pt0->Fill(v0.Pt()/1000.);
                    h_pt1->Fill(v1.Pt()/1000.);

                    double dR12 = v0.DeltaR(v1);
                    h_dr12->Fill(dR12);
                }

                if (i % 500 == 0) {
                    std::cout << "decay length " << dl << " " << dl0 << " ("
                              << x1 << ", " << y1 << ", " << z1 << ")" << '\n';
                }
            }
        }
    }

    new TCanvas("pt");
    h_pt->Draw();

    new TCanvas("eta");
    h_eta->Draw();

    new TCanvas("phi");
    h_phi->Draw();

    (new TCanvas("eta_phi"))->SetLogz(1);
    h_eta_phi->Draw("colz");

    new TCanvas("pt0");
    h_pt0->Draw();

    new TCanvas("pt1");
    h_pt1->Draw();

    new TCanvas("e");
    h_e->Draw();

    new TCanvas("m");
    h_m->Draw();

    new TCanvas("dl");
    h_dl->Draw();
    h_dl->SetXTitle("Decay length [mm]");

    new TCanvas("dl0");
    h_dl0->Draw();
    h_dl0->SetXTitle("Proper lifetime [mm/c]");
    h_dl0->Fit("expo");

    new TCanvas("children");
    h_children->Draw();

    new TCanvas("dr12");
    h_dr12->Draw();
}
