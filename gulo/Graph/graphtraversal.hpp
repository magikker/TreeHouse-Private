#ifndef GULO_GRAPHTRAVERSAL_HPP
#define GULO_GRAPHTRAVERSAL_HPP

#include "forward.hpp"
#include "visitorstate.hpp"

namespace Gulo
{
    namespace Detail
    {
        template <typename DirectedOption>
        class GraphTraversal
        {
        public:

            template <typename Graph, typename Visitor>
            static void visit(Graph & graph, Visitor && visitor)
            {
                typedef typename std::remove_reference<Visitor>::type RawVisitor;
                char * visitorState = (char*)calloc(graph.vertices().size(), sizeof(char)); 
                size_t i = 0;
                while (i < graph.vertices().size()) {
                    while (i < graph.vertices().size() && (graph.vertices()[i] == nullptr || visitorState[i] != VisitorState::Unvisited)) {
                        ++i;
                    }
                    if (!visitor.terminate()) {
                        if (i < graph.vertices().size()) {
                            RawVisitor::Traversal::template visit<DirectedOption>(graph.vertices()[i], visitor, visitorState);
                        }
                    } else {
                        break;
                    }
                }
                free(visitorState);
            }

            template <typename Graph, typename Visitor>
            static void visit(const Graph & graph, Visitor && visitor)
            {
                typedef typename std::remove_reference<Visitor>::type RawVisitor;
                char * visitorState = (char*)calloc(graph.vertices().size(), sizeof(char)); 
                size_t i = 0;
                while (i < graph.vertices().size()) {
                    while (i < graph.vertices().size() && (graph.vertices()[i] == nullptr || visitorState[i] != VisitorState::Unvisited)) {
                        ++i;
                    }
                    if (!visitor.terminate()) {
                        if (i < graph.vertices().size()) {
                            RawVisitor::Traversal::template visit<DirectedOption>(graph.vertices()[i], visitor, visitorState);
                        }
                    } else {
                        break;
                    }
                }
                free(visitorState);
            }
        };
    }
}

#endif // GULO_GRAPHTRAVERSAL_HPP