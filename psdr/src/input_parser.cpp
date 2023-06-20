#include <iostream>
#include <input_parser.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>


namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace std;

int Parser::parse(int argc, char const* argv[]){

    try{

        po::options_description input_options("PSDR");
        input_options.add_options()
                ("help,h", "Help message")
                ("working_dir,w", po::value<string>(), "Working directory.\nAll paths will be treated relative to this directory.")
                ("input_file,i", po::value<string>()->required(), "Input file")
                ("output_file,o", po::value<string>(), "Output file")
                ("source,s", po::value<string>()->default_value("ply"), "Data source:"
                                                                "\n\t-ply"
                                                                "\n\t-npz"
                                                                "\n\t-omvs (an OpenMVS project file)")
                ("groundtruth_file,g", po::value<string>(), "Groundtruth file")
                ("occ", po::value<string>(), "Occupancy file")
                ("transformation_file,t", po::value<string>(), "Transformation file")
                ("crop_file,c", po::value<string>(), "Crop file")
                ("scale", po::value<double>()->default_value(0.0), "Scale mean edge length of Delaunay to this value.")
                ("adt", po::value<double>()->default_value(-1), "Epsilon for adaptive 3DT.")
                ("icomp", po::value<int>()->default_value(1), "Number of connected components of input to keep. -1 = all.")
                ("iclose", po::value<int>()->default_value(1), "Try to close input open meshes with hole filling.")
            ;




        // parse all options to a single options map called vm
        po::store(po::parse_command_line(argc, argv, input_options), vm);

        if (vm.count("help")){
            cout << input_options << "\n";
            return 1;
        }

        // There must be an easy way to handle the relationship between the
        // option "help" and "host"-"port"-"config"
        // Yes, the magic is putting the po::notify after "help" option check
        po::notify(vm);
    }
    catch(std::exception& e){
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    catch(...){
        std::cerr << "Unknown error when parsing command line options!" << "\n";
        return 1;
    }

    return 0;
}
