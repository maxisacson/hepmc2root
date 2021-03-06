#include <iostream>
#include <fstream>
#include <sstream>

#include <unistd.h>

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

void usage(char** argv) {
    std::cout << "Split a single hepmc2-file into several.\n\n";
    std::printf("Usage: %s <input> [base_output] [options]\n", argv[0]);
    std::cout << "Options:\n";
    std::cout << "  -h      Display this message and exit.\n";
    std::cout << "  -n <N>  Process only the first N events.\n";
    std::cout << "  -e <E>  Split input into E events per file.\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage(argv);
        return 1;
    }

    int c;
    int events_per_file = -1;
    int maxevents = -1;
    while ((c = getopt(argc, argv, "hn:e:")) != -1) {
        switch (c) {
            case 'h':
                usage(argv);
                return 0;
                break;
            case 'n':
                {
                    maxevents = atoi(optarg);
                }
                break;
            case 'e':
                {
                    events_per_file = atoi(optarg);
                }
                break;
            default:
                return 2;
                break;
        }
    }

    std::string fn_input = argv[optind];
    std::string fn_output_base = fn_input;
    if (argc >= optind+2) {
        fn_output_base = argv[optind+1];
    }
    std::stringstream fn_output;

    std::ifstream is(fn_input);
    HepMC::GenEvent evt;
    int ievent   = 0;
    int file_ievent   = 0;
    int ifile    = 0;

    HepMC::IO_GenEvent* ascii_io = nullptr;

    std::cout << "Splitting input " << fn_input << " ...\n";
    while (is) {
        if (maxevents >= 0 && ievent >= maxevents) {
            break;
        }

        if (events_per_file > 0 && file_ievent >= events_per_file) {
            file_ievent = 0;
            if (ascii_io != nullptr) {
                delete ascii_io;
                ascii_io = nullptr;
            }
        }

        evt.read(is);

        if (evt.is_valid()) {
            if (ascii_io == nullptr) {
                std::stringstream().swap(fn_output);
                fn_output << fn_output_base << "." << ifile++;
                ascii_io = new HepMC::IO_GenEvent(fn_output.str(), std::ios::out);
            }
            ascii_io->write_event(&evt);
            ++ievent;
            ++file_ievent;
        }
    }
    std::cout << ievent << " events split over " << ifile << " files." << '\n';

    delete ascii_io;

    return 0;
}
