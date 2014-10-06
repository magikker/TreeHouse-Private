#ifndef GULO_GRAPH_WEAKLYCONNECTEDCOMPONENTSOUTDIRECTED_HPP
#define GULO_GRAPH_WEAKLYCONNECTEDCOMPONENTSOUTDIRECTED_HPP

#include <vector>
#include <set>
#include "DepthFirst"

namespace Gulo
{
    namespace Detail
    {
        template <typename VertexPointer>
        struct OutDirectedWeaklyConnectedComponentGatherer : public DepthFirstVisitor
        {
            void startTraversal( VertexPointer v )
            {
                components.push_back({});
                if (componentPos.size() == 0) {
                    componentPos.resize(v->graph()->vertices().size(), -1);
                }
                adjacentComponents.clear();
            }

            void finishTraversal( VertexPointer )
            {
                if (adjacentComponents.size() != 0) {
                    adjacentComponents.insert(components.size()-1);
                    int first = *(std::min_element(adjacentComponents.begin(), adjacentComponents.end()));
                    for (auto i : adjacentComponents) {
                        if (i != first) {
                            for (auto * k : components[i]) {
                                componentPos[k->id()] = first;
                            }
                            components[first].insert(components[i].begin(), components[i].end());
                            components[i].clear();
                        }
                    }
                    components[first].insert(components.back().begin(), components.back().end());
                    components.back().clear();
                } 
            }

            void startVertex(VertexPointer v)
            {
                componentPos[v->id()] = components.size()-1;
                components.back().insert(v);
            }

            void reachVisitedVertex(VertexPointer v)
            {
                if (componentPos[v->id()] != components.size()-1) {
                    adjacentComponents.insert(componentPos[v->id()]);
                }
            }

            std::set<size_t> adjacentComponents;
            std::vector<std::set<VertexPointer> > components;
            std::vector<size_t> componentPos;
        };
    }
}

#endif // GULO_GRAPH_WEAKLYCONNECTEDCOMPONENTSOUTDIRECTED_HPP