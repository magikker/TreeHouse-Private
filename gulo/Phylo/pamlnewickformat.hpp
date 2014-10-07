#ifndef GULO_PAMLNEWICKFORMAT_HPP
#define GULO_PAMLNEWICKFORMAT_HPP

#include <sstream>
#include <queue>
#include <set>
#include "newickformat.hpp"

namespace Gulo
{
    struct PAMLNewickFormat
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
        void setDivergenceDate(BGraph & g, VertexIdentifier n, double date, const std::vector<std::string> & comments)                            // set divergence date
        void setBranchClass(BGraph & graph, EdgeIdentifier edge, int branchclass, const std::vector<std::string> & comments)                    // set branch class        
        void setRoot(MyGraph & g, VertexIdentifier n)                                                                                             // root tree
        void setUnrooted(MyGraph & g)                                                                                                           // unroot tree
        void abort(MyGraph & g)                                                                                                                 // clear memory following parse error
        
        */
        
        template <typename Graph, typename ElementData, typename Mapping>
        struct Processor : public NewickFormat::Processor<Graph, ElementData, Mapping>
        {
            virtual ~Processor()
            {
            }
            
            virtual void processElementData(ElementData & data, Graph & graph, typename Mapping::VertexIdentifier & vertex, typename Mapping::EdgeIdentifier & edge, Mapping & mapping, Newick::ParseResult & result)
            {
                if (!result.hasError()) {
                    if (data.template triggerData<'@'>().size() > 1) {
                        result.setError("Multiple divergence dates specified for the same vertex", data.template triggerData<'@'>()[1].pos());
                        return;
                    } else if (data.template triggerData<'@'>().size() == 1) {
                        mapping.setDivergenceDate(graph, vertex, data.template triggerData<'@'>()[0].value(), data.template triggerData<'@'>()[0].template comments<'['>());
                    }
                    if (data.template triggerData<'#'>().size() > 1) {
                        result.setError("Multiple branch classes specified for the same edge", data.template triggerData<'#'>()[1].pos());
                        return;
                    }
                    if (data.template triggerData<'$'>().size() > 1) {
                        result.setError("Multiple clade classes specified for the same edge", data.template triggerData<'$'>()[1].pos());
                        return;
                    }
                }
            }
            
            virtual void finalize(Graph & graph, std::vector<ElementData> & vertices, std::vector<std::pair<unsigned int, unsigned int> > & edges, Mapping & mapping, Newick::ParseResult & result, std::vector<typename Mapping::VertexIdentifier> & vertexElements, std::vector<typename Mapping::EdgeIdentifier> & edgeElements)
            {
                if (!result.hasError()) {
                    std::map<unsigned int, int> classes;
                    std::map<unsigned int, std::vector<unsigned int> > edges2;
                    for (auto e : edges) {
                        edges2[e.first].push_back(e.second);
                    }
                    for (auto edge : edgeElements) {
                        mapping.setBranchClass(graph, edge, 0, std::vector<std::string>());
                    }
                    auto it = edges.begin();
                    while (it != edges.end()) {
                        if (vertices[it->second].template triggerData<'$'>().size() == 1) {
                            auto e = subtreeEdgeTraversal(it->second, edges2);
                            for (auto pos : e) {
                                mapping.setBranchClass(graph, edgeElements[pos], vertices[it->second].template triggerData<'$'>()[0].value(), vertices[it->second].template triggerData<'$'>()[0].template comments<'['>());
                                classes[pos] = vertices[it->second].template triggerData<'$'>()[0].value();
                            }
                        }
                        ++it;
                    }
                 
                    it = edges.begin();
                    while (!result.hasError() && it != edges.end()) {
                        if (vertices[it->second].template triggerData<'#'>().size() == 1) {
                            mapping.setBranchClass(graph, edgeElements[it->second], vertices[it->second].template triggerData<'#'>()[0].value(), vertices[it->second].template triggerData<'#'>()[0].template comments<'['>());
                            classes[it->second] = vertices[it->second].template triggerData<'#'>()[0].value();                             
                        }
                        ++it;
                    }      
                      
                    std::set<int> finalclasses;
                    for (auto n : classes) {
                        finalclasses.insert(n.second);
                    }
                    classes.erase(0);
                    int z = 1;
                    while (!result.hasError() && z <= finalclasses.size()) {
                        if (finalclasses.find(z) == finalclasses.end()) {
                            std::ostringstream oss;
                            oss << z;
                            result.setError("Expected PAML branch classes to be a consecutive range of integers starting from 0 or 1. Class " + oss.str() + " is missing", 0);                             
                       }
                        ++z;
                    }
                }
            }
           
            std::vector<unsigned int> subtreeEdgeTraversal(unsigned int pos, std::map<unsigned int, std::vector<unsigned int> > & edges)
            {
                std::vector<unsigned int> result;
                std::queue<unsigned int> queue;
                
                result.push_back(pos);
                auto children = edges[pos];
                for (auto c : children) {
                    queue.push(c);
                }
                while (queue.size() > 0) {
                    pos = queue.front();
                    queue.pop();
                    result.push_back(pos);
                    children = edges[pos];
                    for (auto c : children) {
                        queue.push(c);
                    }
                }
                
                return result;
            }
        };
        
        template <typename Mapping>
        using Parser = Newick::Parser<
                                      Mapping,
                                      Newick::VertexLabel<Newick::QuotedOrUnquotedToken>,
                                      Newick::Comment<'[',']', Newick::NestedComment>,
                                      Newick::Trigger<':', Newick::RealToken>,
                                      Newick::Trigger<'@', Newick::RealToken>,
                                      Newick::Trigger<'#', Newick::IntegerToken>,
                                      Newick::Trigger<'$', Newick::IntegerToken>,
                                      Newick::Processor<Processor> 
                                     >;
        
    };
}

#endif // GULO_PAMLNEWICKFORMAT_HPP