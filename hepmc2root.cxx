#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <iterator>

#include <unistd.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenRanges.h"

#include "TTree.h"
#include "TFile.h"

void usage(char** argv) {
    std::printf("Usage: %s <input> [output] [options]\n", argv[0]);
    std::cout << "Options:\n";
    std::cout << "  -h      Display this message and exit.\n";
    std::cout << "  -n <N>  Process only the first N events.\n";
    std::cout << "  -f      Fast/flat mode (no child/vertex/..).\n";
}

template<typename T, typename F>
auto find_idx(const T& v, F&& pred) {
    const auto& pos = std::find_if(std::begin(v), std::end(v), pred);
    return std::distance(std::begin(v), pos);
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

    std::vector<std::vector<int>> children{};
    std::vector<std::vector<int>> parents{};

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

    std::vector<std::vector<int>> vtx_part_in{};
    std::vector<std::vector<int>> vtx_part_out{};
};

struct Output {
    std::unique_ptr<TFile> file{};
    std::unique_ptr<TTree> tree{};
    Event event{};
};

auto fill_particle(const HepMC::GenParticle& p, Event& event) {
    event.pdg_id.push_back(p.pdg_id());
    event.barcode.push_back(p.barcode());
    event.status.push_back(p.status());
    event.is_final_state.push_back((int)(p.status() == 1));

    const auto& momentum = p.momentum();
    event.pt.push_back(momentum.perp());
    event.e.push_back(momentum.e());
    event.m.push_back(momentum.m());
    event.eta.push_back(momentum.eta());
    event.phi.push_back(momentum.phi());

    return event.pdg_id.size() - 1;
}

void fill_vertex(int i, const HepMC::GenVertex& v, Event& event) {
    event.vtx_barcode[i] = v.barcode();

    const auto& position = v.position();
    event.vtx_x[i]       = position.x();
    event.vtx_y[i]       = position.y();
    event.vtx_z[i]       = position.z();
    event.vtx_t[i]       = position.t();
}

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

    event.children.clear();
    event.parents.clear();

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

    event.vtx_part_in.clear();
    event.vtx_part_out.clear();
    event.vtx_part_in_barcode.clear();
    event.vtx_part_out_barcode.clear();
}

void make_output(const std::string& fn_output, Output& output, bool flat) {
    output.file =
            std::unique_ptr<TFile>(TFile::Open(fn_output.c_str(), "RECREATE"));
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

    if (!flat) {
        output.tree->Branch("prod_vtx",          &output.event.prod_vtx);
        output.tree->Branch("decay_vtx",         &output.event.decay_vtx);
        output.tree->Branch("prod_vtx_barcode",  &output.event.prod_vtx_barcode);
        output.tree->Branch("decay_vtx_barcode", &output.event.decay_vtx_barcode);
        output.tree->Branch("children",          &output.event.children);
        output.tree->Branch("parents",           &output.event.parents);
    }

    output.tree->Branch("pt",  &output.event.pt);
    output.tree->Branch("e",   &output.event.e);
    output.tree->Branch("m",   &output.event.m);
    output.tree->Branch("eta", &output.event.eta);
    output.tree->Branch("phi", &output.event.phi);

    if (!flat) {
        output.tree->Branch("vtx_barcode", &output.event.vtx_barcode);
        output.tree->Branch("vtx_x",       &output.event.vtx_x);
        output.tree->Branch("vtx_y",       &output.event.vtx_y);
        output.tree->Branch("vtx_z",       &output.event.vtx_z);
        output.tree->Branch("vtx_t",       &output.event.vtx_t);

        output.tree->Branch(
                "vtx_part_in_barcode",      &output.event.vtx_part_in_barcode);
        output.tree->Branch(
                "vtx_part_out_barcode",     &output.event.vtx_part_out_barcode);
        output.tree->Branch("vtx_part_in",  &output.event.vtx_part_in);
        output.tree->Branch("vtx_part_out", &output.event.vtx_part_out);
    }
}

int process_evt(const HepMC::GenEvent& evt, Event& event, bool flat) {

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
    const auto& vertices  = evt.vertex_range();

    auto n_particles =
            std::distance(std::begin(particles), std::end(particles));
    auto n_vertices = std::distance(std::begin(vertices), std::end(vertices));

    assert(n_particles >= 0);
    assert(n_vertices >= 0);

    event.children.resize(n_particles);
    event.parents.resize(n_particles);
    event.prod_vtx.resize(n_particles);
    event.decay_vtx.resize(n_particles);
    event.prod_vtx_barcode.resize(n_particles);
    event.decay_vtx_barcode.resize(n_particles);

    event.vtx_barcode.resize(n_vertices);
    event.vtx_x.resize(n_vertices);
    event.vtx_y.resize(n_vertices);
    event.vtx_z.resize(n_vertices);
    event.vtx_t.resize(n_vertices);
    event.vtx_part_in_barcode.resize(n_vertices);
    event.vtx_part_out_barcode.resize(n_vertices);
    event.vtx_part_in.resize(n_vertices);
    event.vtx_part_out.resize(n_vertices);

    for (const auto& p : particles) {
        auto ip = fill_particle(*p, event);
        assert(ip >= 0);
        assert(ip < n_particles);

        const auto& prod_vtx = p->production_vertex();
        if (prod_vtx != nullptr && !flat) {
            event.prod_vtx_barcode[ip] = prod_vtx->barcode();

            auto iv = find_idx(vertices, [&prod_vtx](const auto& v) {
                    return prod_vtx->barcode() == v->barcode();
                    });
            assert(iv >= 0);
            assert(iv < n_vertices);

            fill_vertex(iv, *prod_vtx, event);

            event.prod_vtx[ip] = iv;
            event.vtx_part_out_barcode[iv].push_back(p->barcode());
            event.vtx_part_out[iv].push_back(ip);
        } else {
            event.prod_vtx[ip] = -1;
        }

        const auto& end_vtx = p->end_vertex();
        if (end_vtx != nullptr && !flat) {
            event.decay_vtx_barcode[ip] = end_vtx->barcode();

            auto iv = find_idx(vertices, [&end_vtx](const auto& v) {
                    return end_vtx->barcode() == v->barcode();
                    });
            assert(iv >= 0);
            assert(iv < n_vertices);

            fill_vertex(iv, *end_vtx, event);

            event.decay_vtx[ip] = iv;
            event.vtx_part_in_barcode[iv].push_back(p->barcode());
            event.vtx_part_in[iv].push_back(ip);
        } else {
            event.decay_vtx[ip] = -1;
        }
    }

    if (!flat) {
        for (size_t ip = 0; ip < n_particles; ++ip) {
            const auto iv_prod  = event.prod_vtx[ip];
            const auto iv_decay = event.decay_vtx[ip];

            if (iv_decay > -1) {
                for (const auto& child : event.vtx_part_out[iv_decay]) {
                    event.children[ip].push_back(child);
                }
            }

            if (iv_prod > -1) {
                for (const auto& parent : event.vtx_part_in[iv_prod]) {
                    event.parents[ip].push_back(parent);
                }
            }
        }
    }

    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv);
        return 1;
    }

    int c;
    int maxevents = -1;
    bool flat = false;
    while ((c = getopt(argc, argv, "hn:f")) != -1) {
        switch (c) {
        case 'h': {
            usage(argv);
            return 0;
        } break;
        case 'n': {
            maxevents = atoi(optarg);
        } break;
        case 'f': {
            flat = true;
        } break;
        default:
            return 2;
            break;
        }
    }

    std::string fn_input  = argv[optind];
    std::string fn_output = "out.root";
    if (argc >= optind+2) {
        fn_output = argv[optind+1];
    }

    std::cout << "In: " << fn_input << '\n';
    std::cout << "Out: " << fn_output << '\n';
    std::cout << "Processing " << maxevents << " events.\n";

    Output output{};
    make_output(fn_output, output, flat);

    std::ifstream is(fn_input);
    HepMC::GenEvent evt;
    int evt_code = 0;
    int ievent   = 0;
    while (is) {

        if (maxevents >= 0 && ievent >= maxevents) {
            break;
        }

        evt.read(is);

        if (ievent == 0) {
            evt.write_units();
        }

        if (ievent % 500 == 0) {
            std::cout << "ievent " << ievent << '\n';
        }

        if (evt.is_valid()) {
            clear(output.event);
            evt_code = process_evt(evt, output.event, flat);
            output.tree->Fill();
            ++ievent;
        }
    }
    std::cout << ievent << " events processed." << '\n';
    output.file->Write();

    return 0;
}
