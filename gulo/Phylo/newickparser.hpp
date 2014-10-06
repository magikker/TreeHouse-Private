#ifndef GULO_NEWICK_PARSER_HPP
#define GULO_NEWICK_PARSER_HPP

#include <map>
#include <stack>
#include <functional>
#include "../Graph/graphparser.hpp"
#include "newickelementdata.hpp"

namespace Gulo
{
    namespace Newick
    {    
        namespace Detail
        {
            /////////////////
            
            template <char Trigger, typename ElementData, typename Token>
            struct GetTriggerData 
            {
                static void get(const std::string & str, unsigned int & i, ElementData & data, ParseResult & result, const CharacterValidator & validator)
                {
                    data.template activeTriggerData<Trigger>().setValue(std::move(Token::process(str, i, result, validator)));
                }
            };
            
            template <char Trigger, typename ElementData>
            struct GetTriggerData <Trigger, ElementData, VoidToken>
            {
                static void get(const std::string & str, unsigned int & i, ElementData & data, ParseResult & result, const CharacterValidator & validator)
                {
                }
            };
            
            ///////////////////
            
            template <char Trigger, char OpeningTag, char ClosingTag, typename ElementData, typename Token>
            struct GetCommentData 
            {
                static void get(const std::string & str, unsigned int & i, ElementData & data, ParseResult & result, const CharacterValidator & validator)
                {
                    data.template activeTriggerData<Trigger>().template addComment<OpeningTag>(std::move(Token::process(str, i, OpeningTag, ClosingTag, result, validator)));
                }
            };
            
            /////////////////
            
            template <char Trigger, typename ElementData>
            struct AddTriggerData
            {
                static void add(ElementData & data, unsigned int pos)
                {
                    data.template addTrigger<Trigger>(pos);
                }
            };
            
            /////////////////
            
            template <typename ... Tokens>
            struct CommentArray;
            
            template <>
            struct CommentArray<>
            {
                template <typename Map, char Trigger, typename ElementData>
                static void fillFunctionMap(Map & map, CharacterValidator & validator)
                {
                }
            };
            
            template <char OpeningTag, char ClosingTag, typename TokenType, typename ... Tokens>
            struct CommentArray<Comment<OpeningTag, ClosingTag, TokenType>, Tokens...>
            {
                template <typename Map, char Trigger, typename ElementData>
                static void fillFunctionMap(Map & map, CharacterValidator & validator)
                {
                    validator.reserve(OpeningTag);
                    map[std::pair<char,char>(Trigger, OpeningTag)] = &Detail::GetCommentData<Trigger, OpeningTag, ClosingTag, ElementData, TokenType>::get;  
                    CommentArray<Tokens...>::template fillFunctionMap<Map, Trigger, ElementData>(map, validator);
                }
            };
            
            template <typename T, typename ... Tokens>
            struct CommentArray<T, Tokens...>
            {
                template <typename Map, char Trigger, typename ElementData>
                static void fillFunctionMap(Map & map, CharacterValidator & validator)
                {
                    CommentArray<Tokens...>::template fillFunctionMap<Map, Trigger, ElementData>(map, validator);
                }
            };
            
            //////////////////////
            
            template <typename ... Tokens>
            struct TriggerArray;
            
            template <>
            struct TriggerArray<>
            {
                template <typename TriggerFunctionMap, typename TriggerAddFunctionMap,  typename CommentFunctionMap, typename ElementData, typename ... AllTokens>
                static void fillFunctionMaps(TriggerFunctionMap & tmap, TriggerAddFunctionMap & taddmap, CommentFunctionMap & cmap, CharacterValidator & validator)
                {
                    validator.reserve('(');
                    validator.reserve(')');
                    validator.reserve(',');
                    validator.reserve(';');
                }
            };
            
            template <char Tag, typename TokenType, typename ... Tokens>
            struct TriggerArray<Trigger<Tag, TokenType>, Tokens...>
            {
                template <typename TriggerFunctionMap, typename TriggerAddFunctionMap,  typename CommentFunctionMap, typename ElementData, typename ... AllTokens>
                static void fillFunctionMaps(TriggerFunctionMap & tmap, TriggerAddFunctionMap & taddmap, CommentFunctionMap & cmap, CharacterValidator & validator)
                {
                    validator.reserve(Tag);
                    tmap[Tag] = &Detail::GetTriggerData<Tag, ElementData, TokenType>::get;
                    taddmap[Tag] = &Detail::AddTriggerData<Tag, ElementData>::add;
                    CommentArray<AllTokens...>::template fillFunctionMap<CommentFunctionMap, Tag, ElementData>(cmap, validator);
                    TriggerArray<Tokens...>::template fillFunctionMaps<TriggerFunctionMap, TriggerAddFunctionMap, CommentFunctionMap, ElementData, AllTokens...>(tmap, taddmap, cmap, validator);
                }
            };
            
            template <typename TokenType, typename ... Tokens>
            struct TriggerArray<VertexLabel<TokenType>, Tokens...>
            {
                template <typename TriggerFunctionMap, typename TriggerAddFunctionMap,  typename CommentFunctionMap, typename ElementData, typename ... AllTokens>
                static void fillFunctionMaps(TriggerFunctionMap & tmap, TriggerAddFunctionMap & taddmap, CommentFunctionMap & cmap, CharacterValidator & validator)
                {
                    tmap[' '] = &Detail::GetTriggerData<' ', ElementData, TokenType>::get;
                    CommentArray<AllTokens...>::template fillFunctionMap<CommentFunctionMap, ' ', ElementData>(cmap, validator);
                    TriggerArray<Tokens...>::template fillFunctionMaps<TriggerFunctionMap, TriggerAddFunctionMap, CommentFunctionMap, ElementData, AllTokens...>(tmap, taddmap, cmap, validator);
                }
            };
            
            template <typename T, typename ... Tokens>
            struct TriggerArray<T, Tokens...>
            {
                template <typename TriggerFunctionMap, typename TriggerAddFunctionMap, typename CommentFunctionMap, typename ElementData, typename ... AllTokens>
                static void fillFunctionMaps(TriggerFunctionMap & tmap, TriggerAddFunctionMap & taddmap, CommentFunctionMap & cmap, CharacterValidator & validator)
                {
                    TriggerArray<Tokens...>::template fillFunctionMaps<TriggerFunctionMap, TriggerAddFunctionMap, CommentFunctionMap, ElementData, AllTokens...>(tmap, taddmap, cmap, validator);
                }
            };
            
            /////////////////////////////////
        
            // Ensure that the tokens contain no more than one processor item
            
            template <unsigned int I, typename... Tokens>
            struct ProcessorCounter;
            
            template <unsigned int I>
            struct ProcessorCounter<I>
            {
                static const int value = I;
            };
            
            template <unsigned int I, template <typename ...> class TokenType, typename ... Tokens>
            struct ProcessorCounter<I, Processor<TokenType>, Tokens...>
            {
                static const int value = ProcessorCounter<I+1, Tokens...>::value;
            };
            
            template <unsigned int I, typename T, typename ... Tokens>
            struct ProcessorCounter<I, T, Tokens...>
            {
                static const int value = ProcessorCounter<I, Tokens...>::value;
            };
            
            // Get the processor
            
            template <typename Graph, typename ElementData, typename Mapping>
            struct GenericProcessor
            {
                static void process(Graph & graph, 
                            std::vector<ElementData> & vertices, 
                            std::vector<std::pair<unsigned int, unsigned int> > & edges, 
                            Mapping & mapping,
                            ParseResult & result)
                {
                    std::vector<typename Mapping::VertexIdentifier> graphVertices;
                    graphVertices.assign(vertices.size(), typename Mapping::VertexIdentifier());
                    auto it = edges.begin();
                    graphVertices[it->first] = mapping.addVertex(graph, vertices[it->first]);
                    while (it != edges.end()) {
                        graphVertices[it->second] = mapping.addVertex(graph, graphVertices[it->first], vertices[it->second]);  
                        ++it;
                    }
                }
            };
            
            template <typename ... Tokens>
            struct GetProcessor;
            
            template <>
            struct GetProcessor<>
            {
                template <typename Graph, typename ElementData, typename Mapping>
                using type = GenericProcessor<Graph, ElementData, Mapping>;
            };
            
            template <template <typename ...> class TokenType, typename ... Tokens>
            struct GetProcessor<Processor<TokenType>, Tokens...>
            {
                template <typename Graph, typename ElementData, typename Mapping>
                using type = TokenType<Graph, ElementData, Mapping>;
            };
            
            template <typename T, typename ... Tokens>
            struct GetProcessor<T, Tokens...>
            {
                template <typename Graph, typename ElementData, typename Mapping>
                using type = typename GetProcessor<Tokens...>::template type<Graph, ElementData, Mapping>;
            };
        }
        
        //////////////////////
        
        template <typename Mapping, typename... Tokens>
        class Parser
        {            
            template <typename T> using TriggerDatum = Detail::TriggerDatum<T, Tokens..., Trigger<0, VoidToken> >;
            using                       ElementData = Detail::ElementData<Tokens..., Trigger<0, VoidToken> >;
            
            //static_assert(Detail::ProcessorCounter<0, Tokens...>::value <= 1, "Parser requires no more than one Processor template parameter.");
            
            typedef std::map<char, std::function<void(ElementData&, unsigned int)> > TriggerAddFunctionMap;
            typedef std::map<char, std::function<void(const std::string &, unsigned int &, ElementData &, ParseResult &, const CharacterValidator&)> > TriggerFunctionMap;
            typedef std::map<std::pair<char, char>, std::function<void(const std::string &, unsigned int &, ElementData &, ParseResult &, const CharacterValidator&)> > CommentFunctionMap;
              
        public:
            
            using ParseResult = Newick::ParseResult;
            
            Parser()
            {
                  Detail::TriggerArray<Tokens..., Trigger<0, VoidToken> >::template fillFunctionMaps<TriggerFunctionMap, TriggerAddFunctionMap, CommentFunctionMap, ElementData, Tokens..., Trigger<0, VoidToken> >(triggerMap, triggerAddMap, commentMap, validator);
            }
            
            template <typename Graph>
            ParseResult parse(const std::string & newick, Graph & graph)
            {
                typedef typename Detail::GetProcessor<Tokens...>::template type<Graph, ElementData, Mapping> ProcessorType;
                ParseResult result(newick);
                Mapping mapping;
                
                typedef std::pair<unsigned int, unsigned int> uipair;
                
                std::vector<uipair>        edges;
                std::vector<ElementData>   vertices;
                std::stack<uipair>         elems;
                        
                unsigned int i = 0;
                
                vertices.emplace_back();
                triggerAddMap[0](vertices.back(), i);
                elems.push(uipair(0, -1));
                               
                getComments(0, newick, i, vertices[elems.top().first], result);
                unsigned int j = 0;
                while (j < newick.size() && newick[j] != '(' && newick[j] != ')' && newick[j] != ',' && newick[j] != ';') {
                ++j;   
                }
                if (j < newick.size() && newick[j] != '(') {
                    getElementDescriptors(newick, i, vertices[elems.top().first], result);
                }     
                
                bool terminated = false;
                while (!terminated && !result.hasError() && i < newick.size()) {
                    if (std::isspace(newick[i])) {
                        ++i;
                    } else if (newick[i] == '(') {
                        edges.emplace_back(elems.top().first, vertices.size());
                        elems.push(uipair(vertices.size(), edges.size()-1));
                        vertices.emplace_back();
                        triggerAddMap[0](vertices.back(), i);
                        ++i;
                        getComments(0, newick, i, vertices[elems.top().first], result);
                        j = i + 1;
                        while (j < newick.size() && newick[j] != '(' && newick[j] != ')' && newick[j] != ',' && newick[j] != ';') {
                            ++j;   
                        }
                        if (j == newick.size() || newick[j] != '(') {
                            getElementDescriptors(newick, i, vertices[elems.top().first], result); 
                        }
                    } else if (newick[i] == ')') {
                        if (elems.size() <= 1) {
                            result.setError("Descended out of tree (more closing brackets than opening brackets?)", i);
                        } else {
                            elems.pop();
                            ++i;
                            getComments(0, newick, i, vertices[elems.top().first], result);
                            getElementDescriptors(newick, i, vertices[elems.top().first], result);
                        }
                    } else if (newick[i] == ',') {
                        if (elems.size() <= 1) {
                            result.setError("Descended out of tree (more closing brackets than opening brackets?)", i);
                        } else {
                            elems.pop();
                            edges.emplace_back(elems.top().first, vertices.size());
                            elems.push(uipair(vertices.size(), edges.size()-1));
                            vertices.emplace_back();
                            triggerAddMap[0](vertices.back(), i);
                            ++i;                    
                            getComments(0, newick, i, vertices[elems.top().first], result);
                            j = i + 1;
                            while (j < newick.size() && newick[j] != '(' && newick[j] != ')' && newick[j] != ',' && newick[j] != ';') {
                                ++j; 
                            }
                            if (j == newick.size() || newick[j] != '(') {
                                getElementDescriptors(newick, i, vertices[elems.top().first], result);
                            }
                        }
                    } else if (newick[i] == ';') {
                        terminated = true;
                        ++i;
                        j = i;
                        while (i < newick.size() && std::isspace(newick[i])) {
                            ++i;
                        }
                        if (i < newick.size()) {
                            result.setError("Unexpected tokens after semicolon", i);
                        } else {
                            i = j;
                        }
                    } else {
                        result.setError("Unexpected character in newick string", i);
                    }
                }
                
                if (!result.hasError() && !terminated) {
                    result.setError("Expected a terminal semicolon", i);
                } else if (!result.hasError() && elems.size() > 1) {
                    result.setError("Newick string terminated unexpectedly (more opening brackets than closing brackets?)", i);
                }
                
                if (!result.hasError()) {
                    ProcessorType processor;
                    processor.process(graph, vertices, edges, mapping, result);
                }
                if (result.hasError()) {
                    mapping.abort(graph);
                } 
                
                return result;
            }
            
        private:
            
            void getComments(char trigger, const std::string & newick, unsigned int & i, ElementData & elementData, ParseResult & parseResult)
            {
                while (i < newick.size() && std::isspace(newick[i])) {
                    ++i;
                }
                if (!parseResult.hasError() && i != newick.size()) {
                    auto it = commentMap.find(std::pair<char,char>(trigger, newick[i]));
                    while (it != commentMap.end()) {
                        it->second(newick, i, elementData, parseResult, validator);
                        if (parseResult.hasError()) {
                            return;
                        }
                        while (i < newick.size() && std::isspace(newick[i])) {
                            ++i;
                        }
                        it = (i < newick.size()) ? commentMap.find(std::pair<char,char>(trigger,newick[i])) : commentMap.end();
                    }
                }
            }
            
            void getElementDescriptors(const std::string & newick, unsigned int & i, ElementData & elementData, ParseResult & parseResult)
            {
                char lastTrigger = ' ';
                auto it = triggerMap.find(lastTrigger);
                if (it != triggerMap.end()) {
                    it->second(newick, i, elementData, parseResult, validator);
                }
                if (!parseResult.hasError()) {
                    while (i < newick.size() && std::isspace(newick[i])) {
                        ++i;
                    }
                    getComments(lastTrigger, newick, i, elementData, parseResult);
                }
                if (!parseResult.hasError() && i < newick.size()) {
                    auto it = triggerMap.find(newick[i]);
                    while (!parseResult.hasError() && i < newick.size() && it != triggerMap.end()) {
                        lastTrigger = newick[i];
                        triggerAddMap[lastTrigger](elementData, i);
                        ++i;
                        getComments(lastTrigger, newick, i, elementData, parseResult);
                        it->second(newick, i, elementData, parseResult, validator);
                        while (i < newick.size() && std::isspace(newick[i])) {
                            ++i;
                        }
                        getComments(lastTrigger, newick, i, elementData, parseResult);
                        while (i < newick.size() && std::isspace(newick[i])) {
                            ++i;
                        }
                        it = triggerMap.find(newick[i]);
                    }
                }
                if (!parseResult.hasError()) {
                    getComments(lastTrigger, newick, i, elementData, parseResult);
                }
            }
            
            CharacterValidator validator;
            CommentFunctionMap commentMap;
            TriggerFunctionMap triggerMap;
            TriggerAddFunctionMap triggerAddMap;
          
        };
    }
}

#endif // GULO_NEWICK_PARSER_HPP