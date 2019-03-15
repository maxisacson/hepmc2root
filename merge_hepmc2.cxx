#include <iostream>
#include <fstream>
#include <sstream>

#include <unistd.h>

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

void usage(char** argv) {
    std::cout << "Merge several hepmc2-files into a single one.\n\n";
    std::printf("Usage: %s <inputs> [options]\n", argv[0]);
    std::cout << "<inputs>  Space-separated list of inputs to merge.\n";
    std::cout << "Options:\n";
    std::cout << "  -h      Display this message and exit.\n";
    std::cout << "  -o      Output (default: merged.hepmc).\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv);
        return 1;
    }

    int c;
    std::string output = "merged.hepmc";
    while ((c = getopt(argc, argv, "ho:")) != -1) {
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
            default:
                return 2;
                break;
        }
    }

    std::vector<std::string> inputs;
    for (int i = optind; i < argc; ++i) {
        inputs.emplace_back(argv[i]);
    }

    std::cout << "merging:\n";
    for (auto& f : inputs) {
        std::cout << "  " << f;
    }
    std::cout << "\n---> " << output << '\n';

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

            evt.read(is);
            if (evt.is_valid()) {
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
