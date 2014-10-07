#ifndef GULO_DEPTHFIRSTVISITOR_HPP
#define GULO_DEPTHFIRSTVISITOR_HPP

#include "depthfirsttraversal.hpp"

namespace Gulo
{ 
    struct DepthFirst 
    {   
        template <typename Visitor, typename Vertex, typename Edge>
        struct InternalVisitor
        {
            InternalVisitor(Visitor & vis)
            : visitor(vis)
            {
            }
            
            ~InternalVisitor()
            {
            }

            void startTraversal(Vertex * vertex)
            {
                visitor.startTraversal(vertex);
            }
            
            void finishTraversal(Vertex * vertex)
            {
                visitor.finishTraversal(vertex);
            }
            
            void startVertex(Vertex * vertex)
            {
               visitor.startVertex(vertex);
            }
            
            void finishVertex(Vertex * vertex)
            {
                visitor.finishVertex(vertex);
            }
            
            void startEdge(Edge * edge)
            {
                visitor.startEdge(edge);
            }
            
            void finishEdge(Edge * edge)
            {
                visitor.finishEdge(edge);
            }
            
            void reachUnvisitedVertex(Vertex * vertex)
            {
                visitor.reachUnvisitedVertex(vertex);
            }
            
            void reachQueuedVertex(Vertex * vertex)
            {
                visitor.reachQueuedVertex(vertex);
            }
            
            void reachVisitedVertex(Vertex * vertex)
            {
                visitor.reachVisitedVertex(vertex);
            }
            
            void reachBlockedVertex(Vertex * vertex)
            {
                visitor.reachBlockedVertex(vertex);
            }
           
            bool terminate() const
            {
                return visitor.terminate();
            }
            
            Visitor & visitor;
        };
      
        template <typename DirectedOption, typename Visitor, typename Vertex>
        static void visit(Vertex * vertex, Visitor & visitor, char * visitorState = 0)
        {  
            InternalVisitor<Visitor, Vertex, typename Vertex::Graph::Edge> vis(visitor);
            Detail::DepthFirstTraversal::template go<DirectedOption>(vis, vertex, visitorState);
        }
        
        template <typename DirectedOption, typename Visitor, typename Vertex>
        static void visit(const Vertex * vertex, Visitor & visitor, char * visitorState = 0)
        {  
            InternalVisitor<Visitor, const Vertex, const typename Vertex::Graph::Edge> vis(visitor);
            Detail::DepthFirstTraversal::template go<DirectedOption>(vis, vertex, visitorState);
        }
    };
    
    struct DepthFirstVisitor
    {
        typedef DepthFirst Traversal;
        template <typename Vertex> void startTraversal(Vertex * vertex) {}
        template <typename Vertex> void finishTraversal(Vertex * vertex) {}
        template <typename Vertex> void startVertex(Vertex * vertex) {}
        template <typename Vertex> void finishVertex(Vertex * vertex) {}
        template <typename Edge>   void startEdge(Edge * edge) {}
        template <typename Edge>   void finishEdge(Edge * edge) {}
        template <typename Vertex> void reachUnvisitedVertex(Vertex * vertex) {}
        template <typename Vertex> void reachVisitedVertex(Vertex * vertex) {}
        template <typename Vertex> void reachQueuedVertex(Vertex * vertex) {}
        template <typename Vertex> void reachBlockedVertex(Vertex * vertex) {}
        bool terminate() const {return false;}
    };
}

#endif // GULO_DEPTHFIRSTVISITOR_HPP
