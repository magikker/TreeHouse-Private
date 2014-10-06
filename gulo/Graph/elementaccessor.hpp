#ifndef GULO_GRAPH_ELEMENTACCESSOR_HPP
#define GULO_GRAPH_ELEMENTACCESSOR_HPP

namespace Gulo
{
    namespace Detail
    {
        class GraphElementAccessor
        {
        protected:

            template <typename Element>
            static size_t & getElementId(Element * e)
            {
                return e->m_id;
            }
        };
    }
}

#endif // GULO_GRAPH_ELEMENTACCESSOR_HPP