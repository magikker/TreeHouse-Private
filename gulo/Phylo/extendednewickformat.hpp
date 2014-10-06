#ifndef GULO_EXTENDEDNEWICKFORMAT_HPP
#define GULO_EXTENDEDNEWICKFORMAT_HPP

#include "newicktokens.hpp"
#include "newickparser.hpp"
#include "newickprocessorhelper.hpp"

namespace Gulo
{   
    struct ExtendedNewickFormat
    {
        /*
         
        Mapping object should define the following types:
        
        VertexIdentifier
        EdgeIdentifier
        
        and should implement the following functions:
        
        VertexIdentifier addVertex(MyGraph & graph, const std::vector<std::string> & comments)                                                      // add a vertex 
        EdgeIdentifier addEdge(MyGraph & graph, VertexIdentifier from, VertexIdentifier to)                                                         // add an edge
        EdgeIdentifier trailingEdge(MyGraph & graph)                                                                                            // return the identifier used for a trailing edge
        void setLabel(MyGraph & g, VertexIdentifier n, std::string && label, const std::vector<std::string> & comments)                      // set vertex label
        void setLength(MyGraph & g, EdgeIdentifier e, double length, const std::vector<std::string> & comments)                                 // set edge length
        void setHybridVertexType(MyGraph & g, VertexIdentifier n, const std::string & label, const std::vector<std::string> & comments)             // set hybrid vertex type
        void setRoot(MyGraph & g, VertexIdentifier n)                                                                                             // root tree
        void setUnrooted(MyGraph & g)                                                                                                           // unroot tree
        void abort(MyGraph & g)                                                                                                                 // clear memory following parse error
        
        */
        template <typename Graph, typename ElementData, typename Mapping>
        struct Processor
        {
            static void process(Graph & graph, std::vector<ElementData> & vertices, std::vector<std::pair<unsigned int, unsigned int> > & edges, Mapping & mapping, Newick::ParseResult & result)
            {
                //bool rooted = false;
                //bool unrooted = false;
                //typename Mapping::VertexIdentifier root;
                
                std::vector<typename Mapping::VertexIdentifier> vertexElements;
                vertexElements.assign(vertices.size(), typename Mapping::VertexIdentifier());
                std::map<std::string, unsigned int> hybrids;
               
                vertexElements[0] = (mapping.addVertex(graph, vertices[0].template triggerData<0>().template comments<'['>()));
                std::string label = std::move(vertices[0].template triggerData<' '>().value());
                mapping.setLabel(graph, vertexElements[0], std::move(label), vertices[0].template triggerData<' '>().template comments<'['>());
                if (vertices[0].template triggerData<'#'>().size() > 1) {
                    result.setError("More than one hybrid event specified for a single vertex", vertices[0].template triggerData<'#'>()[1].pos());
                } else if (vertices[0].template triggerData<'#'>().size() == 1) {
                    hybrids[vertices[0].template triggerData<'#'>()[0].value()] = 0;
                    mapping.setHybridVertexType(graph, vertexElements[0], vertices[0].template triggerData<'#'>()[0].value(), vertices[0].template triggerData<'#'>()[0].template comments<'['>());
                }
                if (vertices[0].template triggerData<':'>().size() > 1) {
                    result.setError("More than one length specified for a single edge", vertices[0].template triggerData<':'>()[1].pos());
                } else if (vertices[0].template triggerData<':'>().size() == 1) {
                    mapping.setLength(graph, mapping.trailingEdge(graph), vertices[0].template triggerData<':'>()[0].value(), vertices[0].template triggerData<':'>()[0].template comments<'['>());
                }
                
                for (auto & e : edges) {
                    /*unsigned int from = e.first;*/
                    unsigned int to = e.second;
                    if (vertices[e.first].template triggerData<'#'>().size() > 1) {
                        result.setError("More than one hybrid event specified for a single vertex", vertices[e.first].template triggerData<'#'>()[1].pos());
                    } else if (vertices[e.first].template triggerData<'#'>().size() == 1) {
                        auto it = hybrids.find(vertices[e.first].template triggerData<'#'>()[0].value());
                        if (it != hybrids.end()) {
                            e.first = it->second;
                        }
                    }
                    if (vertices[e.second].template triggerData<'#'>().size() > 1) {
                        result.setError("More than one hybrid event specified for a single vertex", vertices[e.second].template triggerData<'#'>()[1].pos());
                    } else if (vertices[e.second].template triggerData<'#'>().size() == 1) {
                        auto it = hybrids.find(vertices[e.second].template triggerData<'#'>()[0].value());
                        if (it == hybrids.end()) {
                            hybrids[vertices[e.second].template triggerData<'#'>()[0].value()] = e.second;
                            vertexElements[e.second] = mapping.addVertex(graph, vertices[e.second].template triggerData<0>().template comments<'['>());
                            mapping.setHybridVertexType(graph, vertexElements[e.second], vertices[e.second].template triggerData<'#'>()[0].value(), vertices[e.second].template triggerData<'#'>()[0].template comments<'['>());
                            std::string label = vertices[e.second].template triggerData<' '>().value();
                            mapping.setLabel(graph, vertexElements[e.second], std::move(label),  vertices[e.second].template triggerData<' '>().template comments<'['>());
                        } else {
                            e.second = it->second;
                        }
                    } else {
                        vertexElements[e.second] = mapping.addVertex(graph, vertices[e.second].template triggerData<0>().template comments<'['>());
                        std::string label = vertices[e.second].template triggerData<' '>().value();
                        mapping.setLabel(graph, vertexElements[e.second], std::move(label), vertices[e.second].template triggerData<' '>().template comments<'['>());
                    }
                    typename Mapping::EdgeIdentifier edge = mapping.addEdge(graph, vertexElements[e.first], vertexElements[e.second]);
                    if (vertices[to].template triggerData<':'>().size() > 1) {
                        result.setError("More than one length specified for a single edge", vertices[e.second].template triggerData<':'>()[1].pos());
                    } else if (vertices[to].template triggerData<':'>().size() == 1) {
                        mapping.setLength(graph, edge, vertices[to].template triggerData<':'>()[0].value(), vertices[to].template triggerData<':'>()[0].template comments<'['>());
                    }
                } 
            
                if (!result.hasError()) {
                    Newick::Detail::FindRootTag<Graph, ElementData, Mapping>::template find<' ', ':', 0>(graph, vertexElements[0], mapping, vertices, vertexElements, edges, result);
                }
            
            }
        };
        
        template <typename Mapping>
        using Parser = Newick::Parser<
                                      Mapping,
                                      Newick::VertexLabel<Newick::QuotedOrUnquotedToken>,
                                      Newick::Comment<'[',']', Newick::NestedComment>,
                                      Newick::Trigger<':', Newick::RealToken>,
                                      Newick::Trigger<'#', Newick::UnquotedToken>,
                                      Newick::Processor<Processor>
                                     >;
    };
}

#endif // GULO_EXTENDEDNEWICKFORMAT_HPP 
