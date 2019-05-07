#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <unistd.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenRanges.h"
#include "HepMC/IO_GenEvent.h"

void usage(char** argv) {
    std::cout << "Remove particles from hepmc2 files.\n\n";
    std::printf("Usage: %s <inputs> [options]\n", argv[0]);
    std::cout << "<inputs>  Space-separated list of inputs to prune.\n";
    std::cout << "Options:\n";
    std::cout << "  -h         Display this message and exit.\n";
    std::cout << "  -o <name>  Output (default: pruned.hepmc).\n";
    std::cout << "  -d <id>    PID to delete from the event record.\n";
    std::cout << "  -k <id>    PID to keep in the event record.\n";
    std::cout << "\n";
    std::cout << "    NOTE: Both -k and -d can be specified multiple times and operate on the abs-value of the PID.\n";
    std::cout << "    NOTE: -k takes precedence over -d.\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv);
        return 1;
    }

    int c;
    std::vector<int> keep_ids = {};
    std::vector<int> remove_ids = {};
    std::string output = "pruned.hepmc";
    while ((c = getopt(argc, argv, "ho:d:k:")) != -1) {
        switch (c) {
            case 'h':
                {
                    usage(argv);
                    return 0;
                }
                break;
            case 'o':
                {
                    output = std::string(optarg);
                }
                break;
            case 'd':
                {
                    remove_ids.push_back(atoi(optarg));
                }
                break;
            case 'k':
                {
                    keep_ids.push_back(atoi(optarg));
                }
                break;
            default:
                return 2;
                break;
        }
    }

    std::vector<std::string> inputs;
    for (int i = optind; i < argc; ++i) {
        inputs.emplace_back(argv[i]);
    }

    std::cout << "pruning:\n";
    for (auto& f : inputs) {
        std::cout << "  " << f;
    }
    std::cout << "\n---> " << output << '\n';

    std::cout << "keeping:";
    for (auto& id : keep_ids) {
        std::cout << "  " << id;
    }
    std::cout << '\n';

    std::cout << "removing:";
    for (auto& id : remove_ids) {
        std::cout << "  " << id;
    }
    std::cout << '\n';

    auto contains = [](std::vector<int> list, int x) {
        return std::find(list.begin(), list.end(), x) != list.end();
    };

    auto* ascii_io = new HepMC::IO_GenEvent(output, std::ios::out);
    int ievent = 0;
    for (auto& input : inputs) {
        std::ifstream is(input);
        HepMC::GenEvent evt;

        std::cout << "<--- " << input << '\n';

        int file_ievent = 0;
        while (is) {

            if (ievent % 1000 == 0) {
                std::cout << ievent << "   [ " << file_ievent << " ]" << '\n';
            }

            std::vector<HepMC::GenParticle*> prune_particles = {};
            std::vector<HepMC::GenVertex*> prune_vertices = {};

            evt.read(is);
            if (evt.is_valid()) {
                for (auto p : evt.particle_range()) {
                    int abs_id = std::abs(p->pdg_id());

                    bool in_remove_list = contains(remove_ids, abs_id);
                    bool in_keep_list = contains(keep_ids, abs_id);
                    bool prune =
                        (in_remove_list && !in_keep_list) ||
                        (keep_ids.size() > 0 && !in_keep_list);

                    if (prune) {
                        prune_particles.push_back(p);
                    }
                }

                while (prune_particles.size() > 0) {
                    auto p = prune_particles.back();
                    auto v_start = p->production_vertex();
                    auto v_end = p->end_vertex();

                    if (v_end != nullptr) {
                        v_end->remove_particle(p);
                    }

                    if (v_start != nullptr) {
                        v_start->remove_particle(p);
                    }

                    delete p;
                    prune_particles.pop_back();
                }

                for (auto v : evt.vertex_range()) {
                    auto n_in = v->particles_in_size();
                    auto n_out = v->particles_out_size();
                    bool prune = (n_in == 0) && (n_out == 0);
                    if (prune) {
                        prune_vertices.push_back(v);
                    }
                }

                while (prune_vertices.size() > 0) {
                    delete prune_vertices.back();
                    prune_vertices.pop_back();
                }

                ascii_io->write_event(&evt);
                ++ievent;
                ++file_ievent;
            }
        }
    }
    std::cout << "processed " << ievent << " events." << '\n';

    delete ascii_io;

    return 0;
}
