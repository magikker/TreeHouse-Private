#ifndef GULO_NEWICK_PROCESSORHELPER_HPP
#define GULO_NEWICK_PROCESSORHELPER_HPP

#include "newickcommon.hpp"

namespace Gulo
{
    namespace Newick
    {
        namespace Detail
        {
            template <char ...>
            struct FindRootTagHelper;
            
            template <>
            struct FindRootTagHelper<>
            {
                template <typename ElementData, typename VertexIdentifier>
                static void find(std::vector<ElementData>& data, std::vector<VertexIdentifier> n, bool  &rooted, bool &unrooted, VertexIdentifier& root, Newick::ParseResult & result)
                {
                }
            };
            
            template <char A, char... B>
            struct FindRootTagHelper<A, B...>
            {
                template <typename ElementData, typename VertexIdentifier>
                static void find(std::vector<ElementData>& vertexElements, std::vector<VertexIdentifier> vertices, bool  &rooted, bool &unrooted, VertexIdentifier& root, Newick::ParseResult & result)
                {
                    findRootFlag(vertexElements, vertices, rooted, unrooted, root, result);
                    if (!result.hasError()) 
                        FindRootTagHelper<B...>::template find(vertexElements, vertices, rooted, unrooted, root, result);
                }
                
                template < typename ElementData, typename VertexIdentifier>
                static void findRootFlag(std::vector<ElementData>& data, std::vector<VertexIdentifier> n, bool  &rooted, bool &unrooted, VertexIdentifier& root, Newick::ParseResult & result)
                {
                    for (unsigned int i = 0; i < data.size(); ++i) {
                        if (!result.hasError())
                            findRootFlag_impl(data[i].template triggerData<A>(), n[i], rooted, unrooted, root, result);
                    }
                }
                
                
                template <typename TriggerData, typename VertexIdentifier>
                static void findRootFlag_impl(TriggerData & d, VertexIdentifier & vertex, bool & rooted, bool & unrooted, VertexIdentifier& root, Newick::ParseResult & result)
                {
                    for (auto n : d.template comments<'['>()) {
                        if (n == "&R" || n == "&r") {
                            if (!rooted && !unrooted) {
                                rooted = true;
                                root = vertex;
                            } else {
                                if (rooted) {
                                    result.setError("More than one root vertex specified", d.pos());
                                } else {
                                    result.setError("Root vertex specified after tree has already been defined as unrooted", d.pos());
                                }
                                return;
                            } 
                        } else if (n == "&U" || n == "&u") {
                            if (!rooted && !unrooted) {
                                unrooted = true;
                            } else {
                                if (unrooted) {
                                    result.setError("Unrooted tree specified more than once", d.pos());
                                } else {
                                    result.setError("Tree defined as unrooted after a root vertex has already been specified", d.pos());
                                }
                                return;
                            }
                        }
                    }
                }
                
                template <typename TriggerData, typename VertexIdentifier>
                static void findRootFlag_impl(std::vector<TriggerData> & d, VertexIdentifier & n, bool & rooted, bool & unrooted, VertexIdentifier& root, Newick::ParseResult & result)
                {
                    for (unsigned int i = 0; i < d.size(); ++i) {
                        if (!result.hasError())
                        findRootFlag_impl(d[i], n, rooted, unrooted, root, result);
                    }
                }
            };
            
             template <typename Graph, typename ElementData, typename Mapping>
             struct FindRootTag
             {
                template <char ... Triggers>
                static void find(Graph & graph, typename Mapping::VertexIdentifier defaultRoot, Mapping & mapping, std::vector<ElementData>& vertices, std::vector<typename Mapping::VertexIdentifier> vertexElements, std::vector<std::pair<unsigned int, unsigned int> > & edges, Newick::ParseResult & result)
                {
                    bool rooted = false;
                    bool unrooted = false;
                    typename Mapping::VertexIdentifier root;
                    FindRootTagHelper<Triggers...>::template find(vertices, vertexElements, rooted, unrooted, root, result);
                    if (!result.hasError()) {
                        if (rooted) {
                            mapping.setRoot(graph, root);
                        } else if (unrooted) {
                            mapping.setUnrooted(graph);
                        } else {
                        unsigned int rootChildCount = 0;
                            for (auto e : edges) {
                                if (e.first == 0) {
                                    ++rootChildCount;
                                }
                            }
                            if (rootChildCount == 2) {
                                mapping.setRoot(graph, defaultRoot);
                            } else {
                                mapping.setUnrooted(graph);
                            }
                        }
                    }
                }
            };
        }
    }
}

#endif // GULO_NEWICK_PROCESSORHELPER_HPP