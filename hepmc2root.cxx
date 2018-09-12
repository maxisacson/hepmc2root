#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <iterator>

#include "HepMC/GenEvent.h"
#include "HepMC/GenRanges.h"

#include "TTree.h"
#include "TFile.h"

void usage(char** argv) {
    std::printf("Usage: %s <input> [output]\n", argv[0]);
}

struct Event {
    int number{0};
    int n_particles{0};
    int n_vertices{0};
    int mpi{0};
    double scale{0.};
    double alphaQCD{0.};
    double alphaQED{0.};
    int id1{0};
    int id2{0};
    int pdf_id1{0};
    int pdf_id2{0};
    double x1{0.};
    double x2{0.};
    double scalePDF{0.};
    double pdf1{0.};
    double pdf2{0.};

    std::vector<double> weights{};

    std::vector<int> pdg_id{};
    std::vector<int> barcode{};
    std::vector<int> status{};
    std::vector<int> is_final_state{};
    std::vector<int> prod_vtx{};
    std::vector<int> decay_vtx{};
    std::vector<int> prod_vtx_barcode{};
    std::vector<int> decay_vtx_barcode{};

    std::vector<double> pt{};
    std::vector<double> e{};
    std::vector<double> m{};
    std::vector<double> eta{};
    std::vector<double> phi{};

    std::vector<int> vtx_barcode{};

    std::vector<double> vtx_x{};
    std::vector<double> vtx_y{};
    std::vector<double> vtx_z{};
    std::vector<double> vtx_t{};

    std::vector<std::vector<int>> vtx_part_in_barcode{};
    std::vector<std::vector<int>> vtx_part_out_barcode{};
};

struct Output {
    std::unique_ptr<TFile> file{};
    std::unique_ptr<TTree> tree{};
    Event event{};
};

void clear(Event& event) {
    event.number      = 0;
    event.n_particles = 0;
    event.n_vertices  = 0;
    event.mpi         = 0;
    event.scale       = 0.;
    event.alphaQCD    = 0.;
    event.alphaQED    = 0.;
    event.id1         = 0;
    event.id2         = 0;
    event.pdf_id1     = 0;
    event.pdf_id2     = 0;
    event.x1          = 0.;
    event.x2          = 0.;
    event.scalePDF    = 0.;
    event.pdf1        = 0.;
    event.pdf2        = 0.;

    event.weights.clear();
    event.pdg_id.clear();
    event.barcode.clear();
    event.status.clear();
    event.is_final_state.clear();
    event.prod_vtx.clear();
    event.decay_vtx.clear();
    event.prod_vtx_barcode.clear();
    event.decay_vtx_barcode.clear();

    event.pt.clear();
    event.e.clear();
    event.m.clear();
    event.eta.clear();
    event.phi.clear();

    event.vtx_barcode.clear();

    event.vtx_x.clear();
    event.vtx_y.clear();
    event.vtx_z.clear();
    event.vtx_t.clear();

    event.vtx_part_in_barcode.clear();
    event.vtx_part_out_barcode.clear();
}

void make_output(const std::string& fn_output, Output& output) {
    output.file = std::unique_ptr<TFile>(TFile::Open(fn_output.c_str(), "RECREATE"));
    output.file->cd();
    output.tree = std::make_unique<TTree>("nominal", "nominal");

    output.tree->Branch("event_number", &output.event.number);
    output.tree->Branch("n_particles",  &output.event.n_particles);
    output.tree->Branch("n_vertices",   &output.event.n_vertices);
    output.tree->Branch("mpi",          &output.event.mpi);
    output.tree->Branch("scale",        &output.event.scale);
    output.tree->Branch("alphaQCD",     &output.event.alphaQCD);
    output.tree->Branch("alphaQED",     &output.event.alphaQED);
    output.tree->Branch("id1",          &output.event.id1);
    output.tree->Branch("id2",          &output.event.id2);
    output.tree->Branch("pdf_id1",      &output.event.pdf_id1);
    output.tree->Branch("pdf_id2",      &output.event.pdf_id2);
    output.tree->Branch("x1",           &output.event.x1);
    output.tree->Branch("x2",           &output.event.x2);
    output.tree->Branch("scalePDF",     &output.event.scalePDF);
    output.tree->Branch("pdf1",         &output.event.pdf1);
    output.tree->Branch("pdf2",         &output.event.pdf2);

    output.tree->Branch("weights",           &output.event.weights);
    output.tree->Branch("pdg_id",            &output.event.pdg_id);
    output.tree->Branch("barcode",           &output.event.barcode);
    output.tree->Branch("status",            &output.event.status);
    output.tree->Branch("is_final_state",    &output.event.is_final_state);
    output.tree->Branch("prod_vtx",          &output.event.prod_vtx);
    output.tree->Branch("decay_vtx",         &output.event.decay_vtx);
    output.tree->Branch("prod_vtx_barcode",  &output.event.prod_vtx_barcode);
    output.tree->Branch("decay_vtx_barcode", &output.event.decay_vtx_barcode);

    output.tree->Branch("pt",  &output.event.pt);
    output.tree->Branch("e",   &output.event.e);
    output.tree->Branch("m",   &output.event.m);
    output.tree->Branch("eta", &output.event.eta);
    output.tree->Branch("phi", &output.event.phi);

    output.tree->Branch("vtx_barcode", &output.event.vtx_barcode);
    output.tree->Branch("vtx_x",       &output.event.vtx_x);
    output.tree->Branch("vtx_y",       &output.event.vtx_y);
    output.tree->Branch("vtx_z",       &output.event.vtx_z);
    output.tree->Branch("vtx_t",       &output.event.vtx_t);

    output.tree->Branch("vtx_part_in_barcode",  &output.event.vtx_part_in_barcode);
    output.tree->Branch("vtx_part_out_barcode", &output.event.vtx_part_out_barcode);
}

int process_evt(const HepMC::GenEvent& evt, Event& event) {
    const auto find_idx = [](const auto& v, const auto& pred) {
        const auto& pos = std::find_if(std::begin(v), std::end(v), pred);
        return std::distance(std::begin(v), pos);
    };

    event.number      = evt.event_number();
    event.n_particles = evt.particles_size();
    event.n_vertices  = evt.vertices_size();
    event.mpi         = evt.mpi();
    event.scale       = evt.event_scale();
    event.alphaQCD    = evt.alphaQCD();
    event.alphaQED    = evt.alphaQED();

    const auto& pdf_info = evt.pdf_info();
    if (pdf_info != nullptr) {
        event.id1      = pdf_info->id1();
        event.id2      = pdf_info->id2();
        event.pdf_id1  = pdf_info->pdf_id1();
        event.pdf_id2  = pdf_info->pdf_id2();
        event.x1       = pdf_info->x1();
        event.x2       = pdf_info->x2();
        event.scalePDF = pdf_info->scalePDF();
        event.pdf1     = pdf_info->pdf1();
        event.pdf2     = pdf_info->pdf2();
    }


    const auto& weights = evt.weights();
    for (const auto& w : weights) {
        event.weights.push_back(w);
    }

    const auto& particles = evt.particle_range();
    const auto& vertices = evt.vertex_range();

    for (const auto& p : particles) {
        event.pdg_id.push_back(p->pdg_id());
        event.barcode.push_back(p->barcode());
        event.status.push_back(p->status());

        event.is_final_state.push_back((int)(p->status() == 1));

        const auto& prod_vtx = p->production_vertex();
        if (prod_vtx != nullptr) {
            event.prod_vtx_barcode.push_back(prod_vtx->barcode());
        } else {
            event.prod_vtx_barcode.push_back(0);
        }

        const auto& end_vtx = p->end_vertex();
        if (end_vtx != nullptr) {
            event.decay_vtx_barcode.push_back(end_vtx->barcode());
        } else {
            event.decay_vtx_barcode.push_back(0);
        }

#if 0
        auto prod_vtx = p->production_vertex();
        if (prod_vtx != nullptr) {
            output.event.prod_vtx.push_back(find_idx(
                    vertices, [bc = prod_vtx->barcode()](const auto& vtx) { return bc == vtx->barcode(); }));
        }

        auto decay_vtx = p->end_vertex();
        if (decay_vtx != nullptr) {
            output.event.decay_vtx.push_back(find_idx(
                    vertices, [bc = decay_vtx->barcode()](const auto& vtx) { return bc == vtx->barcode(); }));
        }
#endif

        const auto& momentum = p->momentum();
        event.pt.push_back(momentum.perp());
        event.e.push_back(momentum.e());
        event.m.push_back(momentum.m());
        event.eta.push_back(momentum.eta());
        event.phi.push_back(momentum.phi());
    }

    for (const auto& v : vertices) {
        event.vtx_barcode.push_back(v->barcode());

        const auto& position = v->position();
        event.vtx_x.push_back(position.x());
        event.vtx_y.push_back(position.y());
        event.vtx_z.push_back(position.z());
        event.vtx_t.push_back(position.t());


        std::vector<int> part_in_barcode{};
        auto p_in = v->particles_in_const_begin();
        for(; p_in != v->particles_in_const_end(); ++p_in) {
            part_in_barcode.push_back((*p_in)->barcode());
        }
        event.vtx_part_in_barcode.push_back(part_in_barcode);

        std::vector<int> part_out_barcode{};
        auto p_out = v->particles_out_const_begin();
        for(; p_out != v->particles_out_const_end(); ++p_out) {
            part_out_barcode.push_back((*p_out)->barcode());
        }
        event.vtx_part_out_barcode.push_back(part_out_barcode);
    }

    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv);
        return 1;
    }

    std::string fn_input = argv[1];
    std::string fn_output = "out.root";
    if (argc >= 3) {
        fn_output = argv[2];
    }

    std::cout << "In: " << fn_input << '\n';
    std::cout << "Out: " << fn_output << '\n';

    Output output{};
    make_output(fn_output, output);

    std::ifstream is(fn_input);
    HepMC::GenEvent evt;
    int evt_code = 0;
    int ievent = 0;
    while (is) {
        evt.read(is);

        if (ievent % 500 == 0) {
            std::cout << "ievent " << ievent << '\n';
        }

        if (evt.is_valid()) {
            clear(output.event);
            evt_code = process_evt(evt, output.event);
            output.tree->Fill();
            ++ievent;
        }


        // if (ievent >= 100) {
        //     break;
        // }
    }
    std::cout << ievent << " events processed." << '\n';
    output.file->Write();

    return 0;
}
