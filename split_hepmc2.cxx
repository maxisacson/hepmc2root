#include <iostream>
#include <fstream>
#include <sstream>

#include <unistd.h>

#include "HepMC/GenEvent.h"
// #include "HepMC/GenRanges.h"
#include "HepMC/IO_GenEvent.h"

int main(int argc, char** argv) {

    int c;
    int events_per_file = -1;
    int maxevents = -1;
    while ((c = getopt(argc, argv, "hn:e:")) != -1) {
        switch (c) {
        case 'h':
            return 0;
            break;
        case 'n': {
            maxevents = atoi(optarg);
        } break;
        case 'e': {
            events_per_file = atoi(optarg);
        } break;
        default:
            return 2;
            break;
        }
    }

    std::string fn_input = argv[optind];
    std::string fn_output_base = "out_events_test.hepmc";
    std::stringstream fn_output;

    std::ifstream is(fn_input);
    HepMC::GenEvent evt;
    int evt_code = 0;
    int ievent   = 0;
    int file_ievent   = 0;
    int ifile    = 0;

    fn_output << fn_output_base << "." << ifile;
    auto* ascii_io = new HepMC::IO_GenEvent(fn_output.str(), std::ios::out);

    while (is) {
        if (maxevents >= 0 && ievent >= maxevents) {
            break;
        }

        if (events_per_file > 0 && file_ievent >= events_per_file) {
            ++ifile;
            file_ievent = 0;
            delete ascii_io;
            std::stringstream().swap(fn_output);
            fn_output << fn_output_base << "." << ifile;
            ascii_io = new HepMC::IO_GenEvent(fn_output.str(), std::ios::out);
        }

        evt.read(is);

        if (evt.is_valid()) {
            auto tmp = &evt;
            ascii_io->write_event(&evt);
            ++ievent;
            ++file_ievent;
        }
    }
    std::cout << ievent << " events processed." << '\n';

    delete ascii_io;

    return 0;
}
