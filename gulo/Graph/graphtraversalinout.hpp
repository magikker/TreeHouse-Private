#ifndef GULO_GRAPHTRAVERSALINOUT_HPP
#define GULO_GRAPHTRAVERSALINOUT_HPP

#include "forward.hpp"
#include "visitorstate.hpp"

namespace Gulo
{
    namespace Detail
    {
        /////////////////////////////
        // Traverse Outgoing Edges //
        /////////////////////////////

        template <template <typename> class DirectedOption>
        struct TraverseOutgoingEdges<DirectedOption<StoreInOutEdges>>
        {
            template <typename Vertex, typename Func>
            static void go(Vertex* v, Func & func, const typename Vertex::Graph::Edge * prev)
            {
                for (auto * edge : v->outgoingEdges()) {
                    func(edge, edge->vertex2());
                }
            }
        };

        template <template <typename> class DirectedOption>
        struct TraverseOutgoingEdges<InOutUndirected<DirectedOption<StoreInOutEdges>>>
        {
            template <typename Vertex, typename Func>
            static void go(Vertex* v, Func & func, const typename Vertex::Graph::Edge * prev)
            {
                for (auto * edge : v->outgoingEdges()) {
                    func(edge, edge->opposite(v));
                }
                for (auto * edge : v->incomingEdges()) {
                    func(edge, edge->opposite(v));
                }
            }
        };

        /////////////////////
        // Graph Traversal //
        /////////////////////

        template <template <typename> class DirectedOption>
        class GraphTraversal<DirectedOption<StoreInOutEdges>>
        {
        public:

            template <typename Graph, typename Visitor>
            static void visit(Graph & graph, Visitor && visitor)
            {
               doVisit<Graph, Visitor, typename Graph::Vertex*, DirectedOption<StoreInOutEdges>>(graph, visitor);
            }

            template <typename Graph, typename Visitor>
            static void visit(const Graph & graph, Visitor && visitor)
            {
                doVisit<const Graph, Visitor, const typename Graph::Vertex*, DirectedOption<StoreInOutEdges>>(graph, visitor);
            }

        protected:

            template <typename Graph, typename Visitor, typename VertexPointer, typename Direction>
            static void doVisit(Graph & graph, Visitor && visitor)
            {
                typedef typename std::remove_reference<Visitor>::type RawVisitor;
                char * visitorState = (char*)calloc(graph.vertices().size(), sizeof(char)); 
                char * precedes = (char*)calloc(graph.vertices().size(), sizeof(char)); 
                size_t i = 0;
                while (i < graph.vertices().size()) {

                    while (i < graph.vertices().size() && (graph.vertices()[i] == nullptr || visitorState[i] != VisitorState::Unvisited)) {
                        ++i;
                    }
                    if (i == graph.vertices().size()) {
                        break;
                    }
                    VertexPointer v = graph.vertices()[i];
                    precedes[v->id()] = 1;
                    VertexPointer w = v;
                    while (v->incomingEdges().size() > 0) {
                        for (auto * k : v->incomingEdges()) {
                            if (precedes[k->vertex1()->id()] == 0) {
                                v = k->vertex1();
                                break;
                            }
                        }
                        if (v == w) {
                            break;
                        } else {
                            w = v;
                            precedes[v->id()] = 1;
                        }
                    }
                    if (!visitor.terminate()) {
                        RawVisitor::Traversal::template visit<Direction>(v, visitor, visitorState);
                    } else {
                        break;
                    }
                }
                free(precedes);
                free(visitorState);           
            }
        };
    }
}

#endif // GULO_GRAPHTRAVERSAL_HPP