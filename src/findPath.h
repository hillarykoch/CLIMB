#ifndef _findPath_H
#define _findPath_H

#include <RcppArmadillo.h>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/dfs.h>
#include <lemon/adaptors.h>
#include "pathenumeration.h"
// prune_enhancer.h used to be included here, 
// there might be some useful functions in here, even though i've eliminated this method

using namespace lemon;
using namespace std;
using namespace Rcpp;
using namespace RcppArmadillo;

void findPath(ListDigraph& gr,
              ListDigraph::Node& src,
              ListDigraph::Node& trg,
              PathEnumeration& enumeration,
              int d,
              int& num_paths,
              vector<int>& all_paths,
              ListDigraph::ArcMap<bool>& filter,
              ListDigraph::Node& curr_node,
              ListDigraph::NodeMap<int>& layer
              )
{
    //std::vector<int> assoc;

    if(num_paths == 0)
    {
        // FIND INITIAL PATH
    	Dfs<ListDigraph> dfs(gr);

         // Run DFS on full graph
    	dfs.run(src, trg);
    	for(ListDigraph::NodeIt n(gr); n != INVALID; ++n) {
    		if(dfs.reached(n)) {
    			enumeration.push_first(n);
    		}
    	}
        num_paths++;

         // Add path to all_paths
        // d+2 = enumeration.len() when the enumeration is full
        for(int i=0; i < d+2; i++) {
            all_paths.push_back(gr.id(enumeration[i]));
        }
        
         // move_curr_node one back in the path
        curr_node = enumeration[d];

         // POP TARGET NODE
        enumeration.pop_last(); // automatically filters between curr_node and target and decrements curr_node's outArcs

        findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter, curr_node, layer);
    } else // IF NUM_PATHS > 0
    {
        // STOPPING RULE
        while(!(enumeration.len() > 0 && enumeration[0] == src && enumeration.outArcs(enumeration[0]) == 0))
        {
            // WHILE THE CURRENT NODE STILL HAS FEASIBLE OUTGOING PATHS
            while(enumeration.outArcs(curr_node) > 0)
            {
                int sz;
                vector<int> temp;
                FilterArcs<ListDigraph> subgraph(gr, filter);
                Dfs<FilterArcs<ListDigraph> > sub_dfs(subgraph);  

                // find new path based on filter
                sub_dfs.run(curr_node, trg);

                // Add to all_paths
                sz = enumeration.len();
                for(int i = 0; i < sz; i++) {
                    all_paths.push_back(gr.id(enumeration[i]));
                }
                for(ListDigraph::NodeIt n(gr); n != INVALID; ++n) {
                    if(sub_dfs.reached(n) && gr.id(n) != gr.id(curr_node)) {
                        temp.push_back(gr.id(n));
                    }
                }
                sort(temp.begin(), temp.end());

                for(auto i = 0; i < temp.size(); i++)
                {
                    all_paths.push_back(temp[i]);

                     for(ListDigraph::NodeIt n(gr); n != INVALID; ++n)
                    {
                        if(gr.id(n) == temp[i])
                        {
                            enumeration.push_last(n);
                        }
                    }
                }
                    
                enumeration.pop_last();
                curr_node = enumeration[d];

                for(ListDigraph::ArcIt a(gr); a != INVALID; ++a) {
                    filter[a] = enumeration.filter(a);
                }

                // delete stuff I don't need anymore
                std::vector<int>().swap(temp);
                findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter, curr_node, layer);
            }

            // ONLY MOVE BACK 1 IF WE ARE NOT ALREADY AT THE SOURCE NODE
            // OTHERWISE, JUST EXIT
            if(curr_node != src)
            {
                // move_curr_node one back in the path
                curr_node = enumeration[layer[curr_node]-1];

                // pop all nodes after curr_node
                for(int i = layer[curr_node]; i < enumeration.len()-1; i++)
                {
                    enumeration.pop_last();
                }

                 // Unfilter all nodes, except the one that curr_node was connected to
                for(ListDigraph::ArcIt a(gr); a != INVALID; ++a)
                {
                    if(!(gr.source(a) == curr_node) && layer[gr.source(a)] >= layer[curr_node])
                    {
                        enumeration.true_filter(a);
                    }
                }

                // Reset outArcs for every node except curr_node
                for(ListDigraph::NodeIt u(gr); u != INVALID; ++u)
                {
                    if(u != curr_node && layer[u] >= layer[curr_node])
                    {
                        enumeration.reset_outArcs(u);
                    }
                }

                 //   UPDATE FILTER IN SUBGRAPH
                for(ListDigraph::ArcIt a(gr); a != INVALID; ++a) {
                    filter[a] = enumeration.filter(a);
                }

                findPath(gr, src, trg, enumeration, d, num_paths, all_paths, filter, curr_node, layer);
            }
        }
    }
}

#endif
