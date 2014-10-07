#ifndef GULO_GRAPHPARSER_HPP
#define GULO_GRAPHPARSER_HPP

namespace Gulo
{
    template <typename Format, typename Graph>
    struct GraphParserMapping;
    
    template <typename Format, typename Graph, typename Mapping = GraphParserMapping<Format, Graph> >
    struct GraphParser
    {
        typedef typename Format::template Parser<Mapping> Parser;
        
        template <typename Data>
        typename Format::template Parser<Mapping>::ParseResult parse(const Data & data, Graph & graph)
        {
            return parser.template parse(data, graph);
        }
        
        Parser parser;
    };
}

#endif // GULO_GRAPHPARSER_HPP