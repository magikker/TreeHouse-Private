#ifndef GULO_DEPTHFIRSTTRAVERSAL_HPP
#define GULO_DEPTHFIRSTTRAVERSAL_HPP

#include <type_traits>
#include "visitorstate.hpp"

// This file provides a generic depth-first traversal algorithm

namespace Gulo
{
    namespace Detail
    {       
        template <typename Directionality>
        struct TraverseOutgoingEdges;

        struct DepthFirstTraversal
        {
        protected:

            template <typename DirectedOption, typename Visitor, typename Vertex>
            static void doDfs(Visitor & visitor, Vertex * vertex, char * visitorState, typename std::conditional<std::is_const<Vertex>::value , const typename Vertex::Graph::Edge *, typename Vertex::Graph::Edge *>::type prev = nullptr)
            {
                typedef typename std::conditional<std::is_const<Vertex>::value, const typename Vertex::Graph::Edge, typename Vertex::Graph::Edge>::type Edge;
                visitor.startVertex(vertex);
                visitorState[vertex->id()] = VisitorState::Queued;
                const auto fn = [&](Edge * e, Vertex * v)->void {
                    visitor.startEdge(e);
                    if (visitorState[v->id()] == VisitorState::Unvisited) {
                        visitor.reachUnvisitedVertex(v);
                        doDfs<DirectedOption>(visitor, v, visitorState, e);
                    } else if (visitorState[v->id()] == VisitorState::Queued) {
                        visitor.reachQueuedVertex(v);
                    } else if (visitorState[v->id()] == VisitorState::Visited) {
                        visitor.reachVisitedVertex(v);
                    } else if (visitorState[v->id()] == VisitorState::Blocked) {
                        visitor.reachBlockedVertex(v);
                    }
                    visitor.finishEdge(e);
                };
                TraverseOutgoingEdges<DirectedOption>::go(vertex, fn, prev);
                visitor.finishVertex(vertex);
                visitorState[vertex->id()] = VisitorState::Visited;
            }
            
        public:

            template <typename DirectedOption, typename Visitor, typename Vertex>
            static void go(Visitor & visitor, Vertex * vertex, char * visitorState = 0)
            {
                bool ownVisitorState = (visitorState == 0);
                if (ownVisitorState) {
                    visitorState = (char*)calloc(vertex->graph()->vertices().size(), sizeof(char)); 
                }
                visitor.startTraversal(vertex);
                visitor.reachUnvisitedVertex(vertex); 
                doDfs<DirectedOption>(visitor, vertex, visitorState);   
                visitor.finishTraversal(vertex);
                if (ownVisitorState) {
                    free(visitorState);
                }
            }
        };
    }
}

#endif // GULO_DEPTHFIRSTTRAVERSAL_HPP