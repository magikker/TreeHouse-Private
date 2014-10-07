#ifndef GULO_GRAPH_DEBUG_HPP
#define GULO_GRAPH_DEBUG_HPP

#include <map>
#include <set>
#include <boost/lexical_cast.hpp>
#include "forward.hpp"

namespace Gulo
{
    namespace Detail
    {
        template <typename Graph>
        void checkGraphElementIds(const Graph & graph)
        {
            size_t i = 0;
            for (auto * e : graph.edges()) {
                if (e != 0 && e->id() != i) {
                    std::cout << "Error:: Edge number " << i << " has id " << e->id() << std::endl;
                }
                ++i;
            }
            i = 0;
            for (auto * v : graph.vertices()) {
                if (v != 0 && v->id() != i) {
                    std::cout << "Error:: Vertex number " << i << " has id " << v->id() << std::endl;
                }
                ++i;
            }
        }

        template <typename Graph>
        typename std::enable_if<std::is_same<typename Graph::Directionality, Undirected>::value>::type checkVertexStorageMatchesGraphData(const Graph & graph)
        {
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertex;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraph;
            
            for (auto * e : graph.edges()) {
                if (e != nullptr) {
                    fromGraph[e->vertex1()].insert(e);
                    fromGraph[e->vertex2()].insert(e);
                }
            }
            for (auto * v : graph.vertices()) {
                if (v != nullptr) {
                    for (auto * e : v->edges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertex[v].insert(e);
                    }
                }

            }
            for (auto * v : graph.vertices()) {
                if(fromVertex[v] != fromGraph[v]) {
                    std::cout << "Error: Edges associated with vertex number " << v->id() << " do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraph[v]) {
                        std::cout << z->id() << " ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertex[v]) {
                        std::cout << z->id() << " ";
                    }
                    std::cout << std::endl << std::endl;
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<Graph::isDirectedGraph() &&  std::is_same<typename Graph::Directionality::EdgeStoragePolicy, StoreOutEdges>::value>::type checkVertexStorageMatchesGraphData(const Graph & graph)
        {
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertex;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraph;
            
            for (auto * e : graph.edges()) {
                if (e != nullptr) {
                    fromGraph[e->vertex1()].insert(e);
                }
            }
            for (auto * v : graph.vertices()) {
                if (v != nullptr) {
                    for (auto * e : v->outgoingEdges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertex[v].insert(e);
                    }
                }
            }
            for (auto * v : graph.vertices()) {
                if(fromVertex[v] != fromGraph[v]) {
                    std::cout << "Error: Edges associated with vertex number " << v->id() << " (" << v->value() << ") do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraph[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertex[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl << std::endl;
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<Graph::isBidirectedGraph() && std::is_same<typename Graph::Directionality::EdgeStoragePolicy, StoreInOutEdges>::value>::type checkVertexStorageMatchesGraphData(const Graph & graph)
        {
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertexIn;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraphIn;
           std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertexOut;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraphOut;

            for (auto * e : graph.edges()) {
                if (e != nullptr) {
                    fromGraphOut[e->vertex1()].insert(e);
                    fromGraphIn[e->vertex2()].insert(e);
                    if (e->isBidirected()) {
                        fromGraphIn[e->vertex1()].insert(e);
                        fromGraphOut[e->vertex2()].insert(e);
                    }
                }
            }
            for (auto * v : graph.vertices()) {
                if (v != nullptr) {
                    for (auto * e : v->outgoingEdges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertexOut[v].insert(e);
                    }

                    for (auto * e : v->incomingEdges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertexIn[v].insert(e);
                    }
                }
            }
            for (auto * v : graph.vertices()) {
                if(fromVertexIn[v] != fromGraphIn[v]) {
                    std::cout << "Error: In Edges associated with vertex number " << v->id() << " (" << v->value() << ") do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraphIn[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertexIn[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl << std::endl;
                }
                if(fromVertexOut[v] != fromGraphOut[v]) {
                    std::cout << "Error: Out Edges associated with vertex number " << v->id() << " do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraphOut[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertexOut[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl << std::endl;
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<Graph::isDirectedGraph() && std::is_same<typename Graph::Directionality::EdgeStoragePolicy, StoreInOutEdges>::value>::type checkVertexStorageMatchesGraphData(const Graph & graph)
        {
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertexIn;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraphIn;
           std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertexOut;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraphOut;

            for (auto * e : graph.edges()) {
                if (e != nullptr) {
                    fromGraphOut[e->vertex1()].insert(e);
                    fromGraphIn[e->vertex2()].insert(e);
                }
            }
            for (auto * v : graph.vertices()) {
                if (v != nullptr) {
                    for (auto * e : v->outgoingEdges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertexOut[v].insert(e);
                    }

                    for (auto * e : v->incomingEdges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertexIn[v].insert(e);
                    }
                }
            }
            for (auto * v : graph.vertices()) {
                if(fromVertexIn[v] != fromGraphIn[v]) {
                    std::cout << "Error: In Edges associated with vertex number " << v->id() << " (" << v->value() << ") do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraphIn[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertexIn[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl << std::endl;
                }
                if(fromVertexOut[v] != fromGraphOut[v]) {
                    std::cout << "Error: Out Edges associated with vertex number " << v->id() << " do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraphOut[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertexOut[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl << std::endl;
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<Graph::isBidirectedGraph() && std::is_same<typename Graph::Directionality::EdgeStoragePolicy, StoreOutEdges>::value>::type checkVertexStorageMatchesGraphData(const Graph & graph)
        {
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromVertex;
            std::map<const typename Graph::Vertex*, std::set<const typename Graph::Edge*> > fromGraph;
            
            for (auto * e : graph.edges()) {
                if (e != nullptr) {
                    fromGraph[e->vertex1()].insert(e);
                    if (e->isBidirected()) {
                        fromGraph[e->vertex2()].insert(e);
                    }
                }
            }
            for (auto * v : graph.vertices()) {
                if (v != nullptr) {
                    for (auto * e : v->outgoingEdges()) {
                        if (e->vertex1() != v && e->vertex2() != v) {
                            std::cout << "Error: Edge number " << e->id() << "is stored within the adjacency list of vertex number " << v->id()
                            << " but the edge things the id of vertex1 is " << e->vertex1()->id() << " and the id of vertex2 is " << e->vertex2()->id() << std::endl << std::endl;
                        }
                        fromVertex[v].insert(e);
                    }
                }
            }
            for (auto * v : graph.vertices()) {
                if(fromVertex[v] != fromGraph[v]) {
                    std::cout << "Error: Edges associated with vertex number " << v->id() << " (" << v->value() << ") do not match those stored in the graph"  << std::endl << std::endl;
                    std::cout << "    According to the graph, this vertex is connected to: ";
                    for (auto * z : fromGraph[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl;
                    std::cout << "    According to the vertex, this vertex is connected to: ";
                    for (auto * z : fromVertex[v]) {
                        std::cout << z->id() << " (" << z->value() << ") ";
                    }
                    std::cout << std::endl << std::endl;
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<Graph::isBidirectedGraph() && std::is_same<typename Graph::Directionality::EdgeStoragePolicy, StoreOutEdges>::value>::type checkBidirectionality(const Graph & graph)
        {
            for (auto * e : graph.edges()) {
                if (e->isBidirected()) {
                    auto & vec1 = e->vertex1()->outgoingEdges();
                    auto & vec2 = e->vertex2()->outgoingEdges();
                    auto it1 = std::find(vec1.begin(), vec1.end(), e);
                    auto it2 = std::find(vec2.begin(), vec2.end(), e);
                    if (it1 == vec1.end()) {
                        std::cout << "Error: bidirectional edge number " << e->id() << " is not stored as an outgoing edge from its vertex1" << std::endl;
                    }
                    if (it2 == vec2.end()) {
                        std::cout << "Error: bidirectional edge number " << e->id() << " is not stored as an outgoing edge from its vertex2" << std::endl;
                    }
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<Graph::isBidirectedGraph() && std::is_same<typename Graph::Directionality::EdgeStoragePolicy, StoreInOutEdges>::value>::type checkBidirectionality(const Graph & graph)
        {
            for (auto * e : graph.edges()) {
                if (e->isBidirected()) {
                    auto & vec1 = e->vertex1()->outgoingEdges();
                    auto & vec2 = e->vertex2()->outgoingEdges();
                    auto it1 = std::find(vec1.begin(), vec1.end(), e);
                    auto it2 = std::find(vec2.begin(), vec2.end(), e);
                    if (it1 == vec1.end()) {
                        std::cout << "Error: bidirectional edge number " << e->id() << " is not stored as an outgoing edge from its vertex1" << std::endl;
                    }
                    if (it2 == vec2.end()) {
                        std::cout << "Error: bidirectional edge number " << e->id() << " is not stored as an outgoing edge from its vertex2" << std::endl;
                    }
                    auto & vec3 = e->vertex1()->incomingEdges();
                    auto & vec4 = e->vertex2()->incomingEdges();
                    auto it3 = std::find(vec3.begin(), vec3.end(), e);
                    auto it4 = std::find(vec4.begin(), vec4.end(), e);
                    if (it3 == vec3.end()) {
                        std::cout << "Error: bidirectional edge number " << e->id() << " is not stored as an incoming edge to its vertex1 [" << e->vertex1()->id() << "]" << std::endl;
                    }
                    if (it4 == vec4.end()) {
                        std::cout << "Error: bidirectional edge number " << e->id() << " is not stored as an incoming edge to its vertex2 [" << e->vertex2()->id() << "]" << std::endl;
                    }
                }
            }
        }

        template <typename Graph>
        typename std::enable_if<!Graph::isBidirectedGraph()>::type checkBidirectionality(const Graph & graph)
        {
        }

        template <typename Graph>
        void checkGraphTopology(const Graph & graph)
        {
            checkGraphElementIds(graph);
            checkVertexStorageMatchesGraphData(graph);
            checkBidirectionality(graph);
        }

        template <typename Graph>
        typename std::enable_if<!Graph::isBidirectedGraph(), std::string>::type dot(const Graph & graph)
        {
            std::string result;
            std::string edge = "--";
            if (graph.isDirectedGraph()) {
                result = "digraph mygraph {\n";
                edge = "->";
            } else {
                result = "graph mygraph {\n";
            }
            for (auto * v : graph.vertices()) {
                if (v != 0) {
                    result += "    " + boost::lexical_cast<std::string>(v->id()) + " [label=\"" + v->label() + "\"];\n";
                }
            }
            for (auto * e : graph.edges()) {
                if (e != 0) {
                    result += "    ";
                    result += boost::lexical_cast<std::string>(e->vertex1()->id()) + " " + edge + " " + boost::lexical_cast<std::string>(e->vertex2()->id()) + ";\n";
                }
            }
            result += "}";
            return result;
        }

        template <typename Graph>
        typename std::enable_if<!Graph::isBidirectedGraph(), std::string>::type dotpos(const Graph & graph)
        {
            std::string result;
            std::string edge = "--";
            if (graph.isDirectedGraph()) {
                result = "digraph mygraph {\n";
                edge = "->";
            } else {
                result = "graph mygraph {\n";
            }
            for (auto * v : graph.vertices()) {
                if (v != 0) {
                    auto x = boost::lexical_cast<std::string>(v->value().first);
                    auto y = boost::lexical_cast<std::string>(v->value().second);
                    result += "    " + boost::lexical_cast<std::string>(v->id()) + " [shape=\"circle\", fixedsize=true,width=\"0.15\", label=\"\", pos=\"" + x + "," + y + "!\"];\n";
                }
            }
            for (auto * e : graph.edges()) {
                if (e != 0) {
                    result += "    ";
                    result += boost::lexical_cast<std::string>(e->vertex1()->id()) + " " + edge + " " + boost::lexical_cast<std::string>(e->vertex2()->id()) + ";\n";
                }
            }
            result += "}";
            return result;
        }

        template <typename Graph>
        typename std::enable_if<Graph::isBidirectedGraph(), std::string>::type dot(const Graph & graph)
        {
            std::string result;
            std::string edge = "->";
                result = "digraph mygraph {\n";

            for (auto * v : graph.vertices()) {
                if (v != 0) {
                    result += "    " + boost::lexical_cast<std::string>(v->id()) + " [label=\"" + v->value() + "\"];\n";
                }
            }
            for (auto * e : graph.edges()) {
                if (e != 0) {
                    result += "    ";
                    result += boost::lexical_cast<std::string>(e->vertex1()->id()) + " " + edge + " " + boost::lexical_cast<std::string>(e->vertex2()->id()) + " ";
                    if (e->isBidirected()) {
                        result += "[dir = \"both\"];";
                    } else {
                        result += ";";
                    }

                    result += ";\n";
                }
            }
            result += "}";
            return result;
        }
    }
}

#endif // GULO_GRAPH_DEBUG_HPP