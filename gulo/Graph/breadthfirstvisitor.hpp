#ifndef GULO_BREADTHFIRSTVISITOR_HPP
#define GULO_BREADTHFIRSTVISITOR_HPP

#include "breadthfirsttraversal.hpp"

namespace Gulo
{    
    struct BreadthFirst 
    {
        template <typename Visitor, typename Vertex, typename Edge>
        struct InternalVisitor
        {
            InternalVisitor(Visitor & vis, size_t vertexCount)
            : queue((std::pair<Vertex*, Edge*>*)malloc(vertexCount * sizeof(std::pair<Vertex*, Edge*>))),
              queueStart(queue),
              queueEnd(queue),
              visitor(vis)
            {
            }
            
            ~InternalVisitor()
            {
                free(queue);
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
                prevEdge = edge;
                visitor.startEdge(edge);
            }
            
            void finishEdge(Edge * edge)
            {
                visitor.finishEdge(edge);
            }
            
            void reachUnvisitedVertex(Vertex * vertex)
            {
                *queueEnd = std::make_pair(vertex, prevEdge);
                ++queueEnd;
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
            
            std::pair<Vertex*, Edge*> & queueFront() const
            {
                return *queueStart;
            }
            
            void queuePop()
            {
                ++queueStart;
            }
            
            bool queueEmpty() const 
            {
                return queueStart == queueEnd;
            }
            
            bool terminate() const
            {
                return visitor.terminate();
            }
            
            std::pair<Vertex*, Edge*> * queue;
            std::pair<Vertex*, Edge*> * queueStart;
            std::pair<Vertex*, Edge*> * queueEnd;
            Visitor & visitor;
            Edge * prevEdge;
        };
      
         template <typename DirectedOption, typename Visitor, typename Vertex>
         static void visit(Vertex * vertex, Visitor & visitor, char * visitorState = 0)
         {  
            InternalVisitor<Visitor, Vertex, typename Vertex::Graph::Edge> vis(visitor, vertex->graph()->vertexCount());
            Detail::BreadthFirstTraversal::template go<DirectedOption>(vis, vertex, visitorState);
         }
         
         template <typename DirectedOption, typename Visitor, typename Vertex>
         static void visit(const Vertex * vertex, Visitor & visitor, char * visitorState = 0)
         {  
            InternalVisitor<Visitor, const Vertex, const typename Vertex::Graph::Edge> vis(visitor, vertex->graph()->vertexCount());
            Detail::BreadthFirstTraversal::template go<DirectedOption>(vis, vertex, visitorState);
         }
    };
    
    struct BreadthFirstVisitor
    {
        typedef BreadthFirst Traversal;
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

#endif // GULO_BREADTHFIRSTVISITOR_HPP
