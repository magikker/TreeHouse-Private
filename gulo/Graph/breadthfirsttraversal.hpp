#ifndef GULO_BREADTHFIRSTTRAVERSAL_HPP
#define GULO_BREADTHFIRSTTRAVERSAL_HPP

#include <type_traits>
#include "visitorstate.hpp"

// This file provides a generic breadth-first traversal algorithm

namespace Gulo
{ 
    namespace Detail
    {       
        struct BreadthFirstTraversal
        {
            template <typename DirectedOption, typename Visitor, typename Vertex>
            static void go(Visitor & visitor, Vertex * vertex, char * visitorState = 0) {
                typedef typename std::conditional<std::is_const<Vertex>::value, const typename Vertex::Graph::Edge, typename Vertex::Graph::Edge>::type Edge;
                bool ownVisitorState = (visitorState == 0);
                if (ownVisitorState) {
                    visitorState = (char*)calloc(vertex->graph()->vertexCount(), sizeof(char)); 
                }
                visitor.startTraversal(vertex);
                visitor.reachUnvisitedVertex(vertex); 
                
                const auto fn = [&](Edge * e, Vertex * v)->void {
                    visitor.startEdge(e);
                    if (visitorState[v->id()] == VisitorState::Unvisited) {
                        visitor.reachUnvisitedVertex(v);
                        visitorState[v->id()] = VisitorState::Queued;
                    } else if (visitorState[v->id()] == VisitorState::Queued) {
                        visitor.reachQueuedVertex(v); 
                    } else if (visitorState[v->id()] == VisitorState::Visited) {
                        visitor.reachVisitedVertex(v);
                    } else if (visitorState[v->id()] == VisitorState::Blocked) {
                        visitor.reachBlockedVertex(v);
                    }
                    visitor.finishEdge(e);
                };
                    
                while (!visitor.queueEmpty() && !visitor.terminate()) {
                    auto pr = visitor.queueFront();
                    visitor.queuePop();
                    visitor.startVertex(pr.first);
                    TraverseOutgoingEdges<DirectedOption>::go(pr.first, fn, pr.second);
                    visitorState[pr.first->id()] = VisitorState::Visited;
                    visitor.finishVertex(pr.first);
                }
                visitor.finishTraversal(vertex);
                if (ownVisitorState) {
                    free(visitorState);
                }
            }
        };
    }
}

#endif // GULO_BREADTHFIRSTTRAVERSAL_HPP
