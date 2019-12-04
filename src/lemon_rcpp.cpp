// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <lemon/list_graph.h>
#include <lemon/dfs.h>
#include <lemon/lgf_reader.h>
#include <lemon/adaptors.h>
#include "findPath.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace lemon;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// Collect output from path enumeration algorithm
// [[Rcpp::export]]
arma::mat cgetPaths(std::string filepath) {
    ListDigraph gr;
    ListDigraph::NodeMap<int> dim(gr);
    ListDigraph::NodeMap<int> label(gr);
    ListDigraph::NodeMap<int> assoc(gr);
    Dfs<ListDigraph> dfs(gr);
    ListDigraph::ArcMap<bool> filter(gr);
    ListDigraph::Node src;
    ListDigraph::Node trg;
    ListDigraph::Node curr_node;
    std::vector<int> all_paths; // instantiate a resizable vector to contain all paths
    std::fstream infile;
    infile.open(filepath);


    // Read in Directed Graph from lgf.txt
    // "attributes" source and target are declared to be nodes and named src and trg
    // dim gives which "layer" the given node lies in
    digraphReader(gr, infile)
        .nodeMap("label", label)
        .nodeMap("dim", dim)
        .nodeMap("assoc", assoc)
        .node("source", src)
        .node("target", trg)
        .run();

    ListDigraph::NodeMap<int> layer(gr);
    for(ListDigraph::NodeIt u(gr); u != INVALID; ++u) {
        layer[u] = dim[u];
    }

    int d = dim[trg]-1;

    // Enumerate all paths using DFS and backtracking
    int num_paths = 0; // when not 0, enter a different section of "findPath" function
    PathEnumeration enumeration(gr, src, trg);

    findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter, curr_node, layer);

    num_paths = (all_paths.size())/(d+2);
    arma::mat out(d+2, num_paths, arma::fill::none);
    for(auto i = 0; i < all_paths.size(); i++) {
        out(i) = all_paths[i];
    }
    return out.t();
}


// Prune paths using all other [(d choose 2) - d + 1] combinations
// [[Rcpp::export]]
arma::uvec crowMatch(arma::mat assoc, arma::mat nonconsec) {
    // do matching
    int n_assoc = assoc.n_rows;
    int n_valid = nonconsec.n_rows;

    // 1 if keep, 0 otherwise
    arma::uvec keepers(n_assoc, arma::fill::zeros);
    for(int i = 0; i < n_assoc; i++) {
        for(int j = 0; j < n_valid; j++) {
            if(assoc(i,0) == nonconsec(j,0) & assoc(i,1) == nonconsec(j,1)) {
                keepers(i) = 1;
                break;
            }
        }
    }
    return keepers;
}

// Get names of an Rcpp::List
// [[Rcpp::export]]
std::vector<std::string> get_list_names(Rcpp::List L) {
    return L.names();
}


// C version of paste0
// [[Rcpp::export]]
std::string cpaste0(std::vector<std::string> str1) {
    std::string pasted;
    for(std::vector<std::string>::iterator it = str1.begin(); it != str1.end(); ++it) {
        pasted += *it;
    }

    return pasted;
}

// perform action just like R's str_split
// [[Rcpp::export]]
arma::mat cstr_split(std::vector<std::string> strings, std::string split) {
    int num_strings = strings.size();
    arma::mat dims(num_strings, 2, arma::fill::none);

    for(auto i = 0; i < num_strings; i++) {
        int num_substr = strings[i].length();
        std::vector<std::string> tmp1;
        std::vector<std::string> tmp2;
        bool found_split = FALSE;

        for(auto j=0; j < num_substr; j++) {
            if(found_split == FALSE) {
                if(strings[i].substr(j,1) != split) {
                    tmp1.push_back(strings[i].substr(j,1));
                } else {
                    found_split = TRUE;
                }
            } else {
                tmp2.push_back(strings[i].substr(j,1));
            }
        }

        // Still need to paste tmp1, tmp2 each into one individual int
        dims(i,0) = std::stoi(cpaste0(tmp1));
        dims(i,1) = std::stoi(cpaste0(tmp2));
    }

    return dims;
}

// Function to pass to .transform()
// [[Rcpp::export]]
int trans_func(double& x) {
    if(x < 0) {
        return(-1);
    } else if(x > 0) {
        return(1);
    } else {
        return(0);
    }
}

// Find which row in a matrix equals some vector
// [[Rcpp::export]]
arma::vec caccept(arma::mat x, arma::colvec y){
    int b = x.n_rows;
    arma::vec out(b, arma::fill::none);

    for(int i = 0; i < b; i++){
        bool vecmatch = arma::approx_equal(x.row(i), y.t(), "absdiff", 0.001);
        if(vecmatch) {
            out(i) = 0;
        } else {
            out(i) = 1;
        }
    }

    return out;
}

// A different and possibly better way to compute the prior_prop
// [[Rcpp::export]]
arma::colvec cget_prior_count(arma::mat red_class,
                              Rcpp::List mus,
                              arma::mat labels,
                              int d,
                              int n,
                              int dist_tol = 0) {
    int n_pairs = mus.size();
    int n_class = red_class.n_rows;
    arma::mat tags(n, n_pairs, arma::fill::zeros);
    arma::colvec m;
    arma::uvec idx;
    arma::uvec tagidx;
    arma::colvec prop_path(n_class, arma::fill::zeros);
    arma::colvec rSums(n);
    arma::uvec outidx;
    arma::colvec onevec;
    arma::rowvec clsvec;

    // Get Combos
    arma::field<arma::mat> comb(n_pairs);
    for(int i = 0; i < n_pairs; i++) {
        comb(i) = Rcpp::as<arma::mat>(mus[i]);
        comb(i).transform([](double val) { return(trans_func(val)); } );
    }

    // Which fits model which pairs of dims
    arma::mat dims = cstr_split(get_list_names(mus), "_");

    // For each possible class
    for(int i = 0; i < n_class; i++) {
        clsvec = red_class.row(i);

        // For each pair of dims
        for(int j = 0; j < n_pairs; j++) {
            m.set_size(comb(j).n_rows);
            idx.set_size(comb(j).n_rows);

            m = caccept(comb(j), { clsvec(dims(j,0)-1), clsvec(dims(j,1)-1) } );
            if(std::all_of(m.begin(), m.end(), [](double x) { return x == 1.0; }))
            {
                continue;
            }

            idx = find(m == 0);
            tagidx = find(labels.col(j) == (idx(0)+1));

            onevec.ones(tagidx.size());
            tags.elem(tagidx+(n*j)) = onevec;
        }

        rSums = sum(tags, 1);
        outidx = find(rSums >= (n_pairs - dist_tol));
        tags.zeros();
        prop_path(i) = (outidx.size());
    }
    return prop_path;
}

// [[Rcpp::export]]
arma::colvec creduce_by_hamming(arma::mat red_class, arma::colvec fullidx, int hamming_tol, int M) {
    int D = red_class.n_cols;
    arma::colvec matchtest(D, arma::fill::ones);
    arma::colvec cur = red_class.row(fullidx(0)).t();
    arma::colvec within_dist(fullidx.size(), arma::fill::zeros);

    for(int i = 0; i < fullidx.size(); i++) {
        matchtest.ones();
        for(int j = 0; j < D; j++) {
            if(cur(j) == red_class(fullidx(i),j)) {
                matchtest(j) = 0;
            }
        }

        if(sum(matchtest) <= hamming_tol) {
            within_dist(i) = 1;
        }
    }
    return within_dist;
}



// Among all reduced classes, get the indices which correspond to the truth
// FOR SIMULATED DATA ONLY
// [[Rcpp::export]]
arma::colvec cget_true_assoc_idx(arma::mat red_class, arma::mat true_assoc) {
    int nassoc = true_assoc.n_rows;
    arma::uvec trueidx;
    arma::colvec m;
    arma::colvec out(nassoc);
    out.ones();
    out = out*(-1);

    int count = 0;
    for(int i = 0; i < nassoc; i++) {
        m = caccept(red_class, true_assoc.row(i).t());
        trueidx = find(m == 0);
        if(trueidx.size() > 0) {
            out(count) = trueidx(0);
            count += 1;
        }
        m.reset();
        trueidx.reset();
    }

    return (out + 1); // adding 1 because indexing is going back to R
}

