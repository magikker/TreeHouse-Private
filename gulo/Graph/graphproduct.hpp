#ifndef GULO_GRAPH_GRAPHPRODUCT_HPP
#define GULO_GRAPH_GRAPHPRODUCT_HPP

namespace Gulo
{
    namespace Detail
    {
        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        DerivedGraph Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::cartesianProduct(const DerivedGraph & other,
                                                                                                       const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier) const
        {
            DerivedGraph g;
            getVector(g.m_vertices).reserve(m_vertices.size() * other.m_vertices.size());
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (size_t j = 0; j < m_vertices.size(); ++j) {
                    if (vertex(j) != nullptr && other.vertex(i) != nullptr) {
                        vertexMultiplier(g.addVertex(), vertex(j), other.vertex(i));
                    } else {
                        getVector(g.m_vertices).push_back(nullptr);
                    }
                }
            }
            size_t offset = 0;
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (auto * e : m_edges) { 
                    if (e != nullptr) {                       
                        auto * ee = g.addEdge(g.vertex(e->vertex1()->id() + offset), g.vertex(e->vertex2()->id() + offset));
                        static_cast<EdgeElementData&>(*ee) = static_cast<const EdgeElementData&>(*e);
                        CopyEdgeExtension<DirectedOption>::template copy<Edge>(ee, e);
                    } else {
                        getVector(g.m_edges).push_back(nullptr);
                    }
                }
                offset += m_vertices.size();
            }
            offset = m_vertices.size();
            for (size_t i = 0; i < m_vertices.size(); ++i) {
                for (auto * e : other.m_edges) {
                    if (e != nullptr) {
                        auto * ee = g.addEdge(g.vertex(i + e->vertex1()->id() * offset), i + g.vertex(e->vertex2()->id() * offset));
                        static_cast<EdgeElementData&>(*ee) = static_cast<const EdgeElementData&>(*e);
                        CopyEdgeExtension<DirectedOption>::template copy<Edge>(ee, e);
                    } else {
                        getVector(g.m_edges).push_back(nullptr);
                    }
                }
            }
            if (elementReindexingIsSuspended() || other.elementReindexingIsSuspended()) {
                g.suspendElementReindexing();
            }
            return g;
        } 

        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        DerivedGraph Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::rootedProduct(const DerivedGraph & other, const Vertex * root,
                                                                                                    const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier) const
        {
            if (root->graph() == &other) {
                return other.rootedProduct(static_cast<const DerivedGraph&>(*this), root, vertexMultiplier);
            }
            DerivedGraph g;
            getVector(g.m_vertices).reserve(m_vertices.size() * other.m_vertices.size());
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (size_t j = 0; j < m_vertices.size(); ++j) {
                    if (vertex(j) != nullptr && other.vertex(i) != nullptr) {
                        vertexMultiplier(g.addVertex(), this->vertex(j), other.vertex(i));
                    } else {
                        getVector(g.m_vertices).push_back(nullptr);
                    }
                }
            }
            size_t offset = 0;
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (auto * e : this->m_edges) {            
                    if (e != nullptr) {            
                        auto * ee = g.addEdge(g.vertex(e->vertex1()->id() + offset), g.vertex(e->vertex2()->id() + offset));
                        static_cast<EdgeElementData&>(*ee) = static_cast<const EdgeElementData&>(*e);
                        CopyEdgeExtension<DirectedOption>::template copy<Edge>(ee, e);
                    } else {
                        getVector(g.m_edges).push_back(nullptr);
                    }
                }
                offset += this->m_vertices.size();
            }
            offset = this->m_vertices.size();
            for (auto * e : other.m_edges) {
                if (e != nullptr) {
                    auto * ee = g.addEdge(g.vertex(root->id() + e->vertex1()->id() * offset), root->id() + g.vertex(e->vertex2()->id() * offset));
                    static_cast<EdgeElementData&>(*ee) = static_cast<const EdgeElementData&>(*e);
                    CopyEdgeExtension<DirectedOption>::template copy<Edge>(ee, e);
                } else {
                    getVector(g.m_edges).push_back(nullptr);
                }
            }
            if (elementReindexingIsSuspended() || other.elementReindexingIsSuspended()) {
                g.suspendElementReindexing();
            }
            return g;
        }

        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        DerivedGraph Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::lexicographicalProduct(const DerivedGraph & other,
                                                                                                             const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier) const
        {
            DerivedGraph g;
            getVector(g.m_vertices).reserve(m_vertices.size() * other.m_vertices.size());
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (size_t j = 0; j < m_vertices.size(); ++j) {
                    if (vertex(j) != nullptr && other.vertex(i) != nullptr) {
                        vertexMultiplier(g.addVertex(), vertex(j), other.vertex(i));
                    } else {
                        getVector(g.m_vertices).push_back(nullptr);
                    }
                }
            }
            size_t offset = 0;
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (auto * e : m_edges) { 
                    size_t offset2 = 0;
                    for (size_t j = 0; j < other.m_vertices.size(); ++j) {
                        if (e != nullptr) {                       
                            auto * ee = g.addEdge(g.vertex(e->vertex1()->id() + offset), g.vertex(e->vertex2()->id() + offset2));
                            static_cast<EdgeElementData&>(*ee) = static_cast<const EdgeElementData&>(*e);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(ee, e);
                        } else {
                            getVector(g.m_edges).push_back(nullptr);
                        }
                        offset2 += m_vertices.size();
                    }
                }
                offset += m_vertices.size();
            }
            offset = m_vertices.size();
            for (size_t i = 0; i < m_vertices.size(); ++i) {
                for (auto * e : other.m_edges) {
                    if (e != nullptr) {
                        auto * ee = g.addEdge(g.vertex(i + e->vertex1()->id() * offset), i + g.vertex(e->vertex2()->id() * offset));
                        static_cast<EdgeElementData&>(*ee) = static_cast<const EdgeElementData&>(*e);
                        CopyEdgeExtension<DirectedOption>::template copy<Edge>(ee, e);
                    } else {
                        getVector(g.m_edges).push_back(nullptr);
                    }
                }
            }
            if (elementReindexingIsSuspended() || other.elementReindexingIsSuspended()) {
                g.suspendElementReindexing();
            }
            return g;
        }  

        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        DerivedGraph Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::strongProduct(const DerivedGraph & other,
                                                                                                    const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier,
                                                                                                    const std::function<void(DerivedEdge * target, const DerivedEdge * e1, const DerivedEdge * e2)> & edgeMultiplier) const
        {
            DerivedGraph g = cartesianProduct(other, vertexMultiplier);
            for (auto * e1 : m_edges) {
                if (e1 != nullptr) {
                    for (auto * e2 : other.m_edges) {
                        if (e2 != nullptr) {
                            size_t v1 = e2->vertex1()->id() * m_vertices.size() + e1->vertex1()->id();
                            size_t v2 = e2->vertex2()->id() * m_vertices.size() + e1->vertex2()->id();
                            auto * e = g.addEdge(g.vertex(v1), g.vertex(v2));
                            static_cast<EdgeElementData&>(*e) = static_cast<const EdgeElementData&>(*e1);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(e, e1);
                            edgeMultiplier(e, e1, e2);
                            v1 = e2->vertex1()->id() * m_vertices.size() + e1->vertex2()->id();
                            v2 = e2->vertex2()->id() * m_vertices.size() + e1->vertex1()->id();
                            e = g.addEdge(g.vertex(v1), g.vertex(v2));
                            static_cast<EdgeElementData&>(*e) = static_cast<const EdgeElementData&>(*e1);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(e, e1);
                            edgeMultiplier(e, e1, e2);
                        }
                    }
                }
            }
            return g;
        } 

        template <typename DerivedGraph, typename DerivedVertex, typename DerivedEdge, typename DirectedOption>
        DerivedGraph Graph<DerivedGraph, DerivedVertex, DerivedEdge, DirectedOption>::tensorProduct(const DerivedGraph & other,
                                                                                                    const std::function<void(DerivedVertex * target, const DerivedVertex * v1, const DerivedVertex * v2)> & vertexMultiplier,
                                                                                                    const std::function<void(DerivedEdge * target, const DerivedEdge * e1, const DerivedEdge * e2)> & edgeMultiplier) const
        {
            DerivedGraph g;
            getVector(g.m_vertices).reserve(m_vertices.size() * other.m_vertices.size());
            for (size_t i = 0; i < other.m_vertices.size(); ++i) {
                for (size_t j = 0; j < m_vertices.size(); ++j) {
                    if (vertex(j) != nullptr && other.vertex(i) != nullptr) {
                        vertexMultiplier(g.addVertex(), vertex(j), other.vertex(i));
                    } else {
                        getVector(g.m_vertices).push_back(nullptr);
                    }
                }
            }
            for (auto * e1 : m_edges) {
                if (e1 != nullptr) {
                    for (auto * e2 : other.m_edges) {
                        if (e2 != nullptr) {
                            size_t v1 = e2->vertex1()->id() * m_vertices.size() + e1->vertex1()->id();
                            size_t v2 = e2->vertex2()->id() * m_vertices.size() + e1->vertex2()->id();
                            auto * e = g.addEdge(g.vertex(v1), g.vertex(v2));
                            static_cast<EdgeElementData&>(*e) = static_cast<const EdgeElementData&>(*e1);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(e, e1);
                            edgeMultiplier(e, e1, e2);
                            v1 = e2->vertex1()->id() * m_vertices.size() + e1->vertex2()->id();
                            v2 = e2->vertex2()->id() * m_vertices.size() + e1->vertex1()->id();
                            e = g.addEdge(g.vertex(v1), g.vertex(v2));
                            static_cast<EdgeElementData&>(*e) = static_cast<const EdgeElementData&>(*e1);
                            CopyEdgeExtension<DirectedOption>::template copy<Edge>(e, e1);
                            edgeMultiplier(e, e1, e2);
                        }
                    }
                }
            }
            if (elementReindexingIsSuspended() || other.elementReindexingIsSuspended()) {
                g.suspendElementReindexing();
            }
            return g;
        }
    }
}


#endif // GULO_GRAPH_GRAPHPRODUCT_HPP