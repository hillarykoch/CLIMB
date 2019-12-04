#ifndef _PATHENUMERATION_H
#define _PATHENUMERATION_H

#include "llist.h"
#include <lemon/list_graph.h>

class PathEnumeration {
    private:
        const lemon::ListDigraph& _gr;
        lemon::ListDigraph::Node _src;
        lemon::ListDigraph::Node _trg;
        lemon::ListDigraph::NodeMap<int> _deg;
        llist _pathlist;
    	lemon::ListDigraph::ArcMap<bool> _filter;

    public:
        PathEnumeration(const lemon::ListDigraph& gr, lemon::ListDigraph::Node& src, lemon::ListDigraph::Node& trg) :_gr(gr), _deg(gr), _filter(gr) {
            for (lemon::ListDigraph::NodeIt u(gr); u != lemon::INVALID; ++u) {
                _deg[u] = countOutArcs(gr, u);
            }
            for (lemon::ListDigraph::ArcIt a(gr); a != lemon::INVALID; ++a){
            	_filter[a] = true; // initialize every path to be valid at start	
            }
      	}

        int outArcs(const lemon::ListDigraph::Node& u) const {
                return _deg[u];
        }
        
        bool filter(const lemon::ListDigraph::Arc& u) const {
            return _filter[u];
        }
        
        const int len() {
            return _pathlist.len();
        }
        
        const lemon::ListDigraph::Node& operator[](int index) const {
            return _pathlist[index];
        } 
        
        void decrement_arcs(lemon::ListDigraph::Node& u) {
            _deg[u]--;
        }

        void true_filter(lemon::ListDigraph::Arc& a) {
            _filter[a] = true;
        }
        
        void reset_outArcs(const lemon::ListDigraph::Node& u) {
            _deg[u] = lemon::countOutArcs(_gr, u);      
        }
          
        void push_first(const lemon::ListDigraph::Node& u) {
            _pathlist.insert_first(u);
        }
        
        void push_last(const lemon::ListDigraph::Node& u) {
            _pathlist.insert_last(u);
        }
        
        void pop_last() {
            // have not added something to handle when sz < 2
            // though this should never happen here
            int sz = _pathlist.len();
            const lemon::ListDigraph::Node u = _pathlist[sz-1]; // The tail
            const lemon::ListDigraph::Node v = _pathlist[sz-2]; // The predecessor to the tail
        
            for (lemon::ListDigraph::OutArcIt a(_gr, v); a != lemon::INVALID; ++a) {
                if(_gr.source(a) == v && _gr.target(a) == u) {
                    _filter[a] = false; // the edge is hidden in the subgraph
                    _deg[_gr.source(a)]--;
                }
            }
            _pathlist.remove_last();
        }
        
        void pop_first() {
            // have not added something to handle when sz < 2
            // though this should never happen here
            _pathlist.remove_first();
        }
};

#endif