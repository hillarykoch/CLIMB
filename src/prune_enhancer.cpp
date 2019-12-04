// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <lemon/list_graph.h>
#include <lemon/dfs.h>
#include <lemon/lgf_reader.h>
#include <lemon/adaptors.h>
#include <iostream>
#include <fstream>
#include <string>
#include <array>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace lemon;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// test reading in a table
// [[Rcpp::export]]
arma::mat cassociate(arma::mat& paths,
                            std::string filepath,
                            int len_filt_h) {
    // path is a matrix of putative paths
    // filepath is location of node file
    // filt_h is only needed for its total length

    ifstream stream (filepath);
    std::string line;
    //int sz = path.size()-2;
    int sz = paths.n_cols - 2;
    int num_paths = paths.n_rows;
    int x, y, z;
    std::vector<int> c(len_filt_h + 2);
    //std::vector<int> assoc(sz);
    arma::mat assoc(num_paths, sz);

    unsigned int count = 0;
    while(std::getline(stream, line))
    {
        ++count;
        if (count > (len_filt_h + 3)) { break; }
        if(count < 2) { continue; }

        stream>>x>>y>>z;
        c[count-2] = z;
    }

    // this range is because we dont need dummy nodes "src" and "trg"
    for(auto i = 0; i < num_paths; i++)
    {
        for(auto j = 0; j < sz; j++)
        {
            assoc(i,j) = c[paths(i,j+1)];
        }
    }

    stream.close();

    return assoc;
}