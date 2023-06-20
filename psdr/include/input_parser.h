#pragma once

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/uuid.hpp>

using namespace std;
namespace po = boost::program_options;


class Parser{

public:

    po::variables_map vm; // this is a std::map which will be populate with all options available for the current mode
    int description_width = 160;

    int parse(int argc, char const *argv[]); // parse the input args and populate vm map


};

