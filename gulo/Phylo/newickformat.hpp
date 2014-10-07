#ifndef GULO_NEWICKFORMAT_HPP
#define GULO_NEWICKFORMAT_HPP

#include "newicktokens.hpp"
#include "newickparser.hpp"
#include "newickprocessorhelper.hpp"

namespace Gulo
{
    struct NewickFormat
    {
        /*
         
        Mapping object should define the following types:
        
        VertexIdentifier
        EdgeIdentifier
        
        and should implement the following functions:
        
        VertexIdentifier addVertex(MyGraph & graph, const std::vector<std::string> & comments)                                                      // add a vertex    
        EdgeIdentifier addEdge(MyGraph & graph, VertexIdentifier from, VertexIdentifier & to)                                                       // add an edge
        EdgeIdentifier trailingEdge(MyGraph & graph)                                                                                            // return the identifier used for a trailing edge
        void setLabel(MyGraph & g, VertexIdentifier n, std::string && label, const std::vector<std::string> & comments)                      // set vertex label
        void setLength(MyGraph & g, EdgeIdentifier e, double length, const std::vector<std::string> & comments)                                 // set edge length
        void setRoot(MyGraph & g, VertexIdentifier n)                                                                                             // root tree
        void setUnrooted(MyGraph & g)                                                                                                           // unroot tree
        void abort(MyGraph & g)                                                                                                                 // clear memory following parse error
        
        */
        
        template <typename Graph, typename ElementData, typename Mapping>
        struct Processor
        {
            virtual ~Processor()
            {
            }
            
            virtual void initialize(Graph & graph, std::vector<ElementData> & vertices, std::vector<std::pair<unsigned int, unsigned int> > & edges, Mapping & mapping, Newick::ParseResult & result)
            {
            }
            
            virtual void finalize(Graph & graph, std::vector<ElementData> & vertices, std::vector<std::pair<unsigned int, unsigned int> > & edges, Mapping & mapping, Newick::ParseResult & result, std::vector<typename Mapping::VertexIdentifier> & vertexElements, std::vector<typename Mapping::EdgeIdentifier> & edgeElements)
            {
            }
            
            virtual void processElementData(ElementData & data, Graph & graph, typename Mapping::VertexIdentifier & vertex, typename Mapping::EdgeIdentifier & edge, Mapping & mapping, Newick::ParseResult & result)
            {
            }
            
            void process(Graph & graph, std::vector<ElementData> & vertices, std::vector<std::pair<unsigned int, unsigned int> > & edges, Mapping & mapping, Newick::ParseResult & result)
            {
                initialize(graph, vertices, edges, mapping, result);
                
                std::vector<typename Mapping::VertexIdentifier> vertexElements;
                std::vector<typename Mapping::EdgeIdentifier> edgeElements;
                vertexElements.assign(vertices.size(), typename Mapping::VertexIdentifier());
                edgeElements.assign(vertices.size(), typename Mapping::EdgeIdentifier());
               
                auto it = edges.begin();             
                vertexElements[it->first] = mapping.addVertex(graph, vertices[it->first].template triggerData<0>().template comments<'['>());
                edgeElements[it->first] = mapping.trailingEdge(graph);
                // ensure that edge created directly after node and before any data is added, so that tree-restricted objects can create objects before assigning data
                std::string label = vertices[it->first].template triggerData<' '>().value();
                if (label.size() > 0) {
                    mapping.setLabel(graph, vertexElements[it->first], std::move(label), vertices[it->first].template triggerData<' '>().template comments<'['>());
                }
                if (vertices[it->first].template triggerData<':'>().size() > 1) {
                    result.setError("Multiple lengths specified for the same edge", vertices[it->first].template triggerData<':'>()[1].pos());
                    return;
                } else if (vertices[it->first].template triggerData<':'>().size() == 1) {
                    mapping.setLength(graph, edgeElements[it->first], vertices[it->first].template triggerData<':'>()[0].value(), vertices[it->first].template triggerData<':'>()[0].template comments<'['>());
                }
                processElementData(vertices[it->first], graph, vertexElements[it->first], edgeElements[it->first], mapping, result);
   
                while (!result.hasError() && it != edges.end()) {
                    vertexElements[it->second] = mapping.addVertex(graph, vertices[it->second].template triggerData<0>().template comments<'['>());
                    edgeElements[it->second] = mapping.addEdge(graph, vertexElements[it->first], vertexElements[it->second]);
                    std::string label = vertices[it->second].template triggerData<' '>().value();
                    if (label.size() > 0) {
                        mapping.setLabel(graph, vertexElements[it->second], std::move(label), vertices[it->second].template triggerData<' '>().template comments<'['>());
                    }
                    if (vertices[it->second].template triggerData<':'>().size() > 1) {
                        result.setError("Multiple lengths specified for the same edge", vertices[it->second].template triggerData<':'>()[1].pos());
                        return;
                    } else if (vertices[it->second].template triggerData<':'>().size() == 1) {
                        mapping.setLength(graph, edgeElements[it->second], vertices[it->second].template triggerData<':'>()[0].value(), vertices[it->second].template triggerData<':'>()[0].template comments<'['>());
                    }
                    processElementData(vertices[it->second], graph, vertexElements[it->second], edgeElements[it->second], mapping, result);
                    ++it;
                }
                
                if (!result.hasError()) {
                    Newick::Detail::FindRootTag<Graph, ElementData, Mapping>::template find<' ', ':', 0>(graph, vertexElements[0], mapping, vertices, vertexElements, edges, result);
                }
                
                finalize(graph, vertices, edges, mapping, result, vertexElements, edgeElements);
            }
            
        };
        
        template <typename Mapping>
        using Parser = Newick::Parser<
                                      Mapping,
                                      Newick::VertexLabel<Newick::QuotedOrUnquotedToken>,
                                      Newick::Comment<'[',']', Newick::NestedComment>,
                                      Newick::Trigger<':', Newick::RealToken>,
                                      Newick::Processor<Processor>
                                     >;
    };
}

#endif // GULO_NEWICKFORMAT_HPP