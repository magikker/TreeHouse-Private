#ifndef GULO_NEWICK_ELEMENTDATA_HPP
#define GULO_NEWICK_ELEMENTDATA_HPP

#include <tuple>
#include <vector>
#include <type_traits>

#include "newickcommon.hpp"

namespace Gulo
{
    namespace Newick
    {        
        namespace Detail
        {
            template <typename T, unsigned int I = 0>
            struct TriggerCommentData;
            
            template <unsigned int I>
            struct TriggerCommentData<std::tuple<>, I>
            {
                typedef std::tuple<> tuple_type;
                template <char C> using tuple_element_type = void;
                template <char C> using token_type = void;
                template <char C> using index = Index<I>;
            };
            
            template <char OpeningTag, char ClosingTag, typename TokenType, typename ... Tokens, unsigned int I>
            struct TriggerCommentData<std::tuple<Comment<OpeningTag, ClosingTag, TokenType>, Tokens...>, I>
            {
                template <char C> using token_type = typename std::conditional<C == OpeningTag, typename TokenType::type, typename TriggerCommentData<std::tuple<Tokens...>, I+1>::template token_type<C> >::type;                
                template <char C> using tuple_element_type = typename std::conditional<C == OpeningTag, std::vector<typename TokenType::type>, typename TriggerCommentData<std::tuple<Tokens...>, I+1>::template tuple_element_type<C> >::type;                
                using tuple_type = typename Prepend<tuple_element_type<OpeningTag>, typename TriggerCommentData<std::tuple<Tokens...>, I+1>::tuple_type>::type;
                template <char C> using index = typename std::conditional<C == OpeningTag, Index<I>, typename TriggerCommentData<std::tuple<Tokens...>, I+1>::template index<C> >::type;
            };
            
            template <typename T, typename ... Tokens, unsigned int I>
            struct TriggerCommentData<std::tuple<T, Tokens...>, I>
            {
                template <char C> using token_type = typename TriggerCommentData<std::tuple<Tokens...>, I >::template token_type<C>;
                template <char C> using tuple_element_type = typename TriggerCommentData<std::tuple<Tokens...>, I >::template tuple_element_type<C>;
                using tuple_type = typename TriggerCommentData<std::tuple<Tokens...>, I >::tuple_type;
                template <char C> using index = typename TriggerCommentData<std::tuple<Tokens...>, I>::template index<C>;
            };
            
            ///////////////////////
            
            template <typename, typename ...> class TriggerDatum;
            
            ///////////////////////
            
            template <typename T, unsigned int I, typename ... AllTokens>
            struct TriggerData;
            
            template <unsigned int I, typename ... AllTokens>
            struct TriggerData<std::tuple<>, I, AllTokens...>
            {
                template <char C> using token_type = void;
                template <char C> using datum_type = void;
                template <char C> using tuple_element_type = void;
                typedef std::tuple<> tuple_type;
                template <char C> using index = Index<I>;
            };
            
            template <char Tag, typename TokenType, typename ... Tokens, unsigned int I, typename ... AllTokens>
            struct TriggerData<std::tuple<Trigger<Tag, TokenType>, Tokens...>, I, AllTokens... >
            {
                template <char C> using token_type = typename std::conditional<C == Tag, typename TokenType::type, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens...>::template token_type<C> >::type;
                template <char C> using datum_type = typename std::conditional<C == Tag, TriggerDatum<typename TokenType::type, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens...>::template datum_type<C> >::type;
                template <char C> using tuple_element_type = typename std::conditional<C == Tag, std::vector<TriggerDatum<typename TokenType::type, AllTokens...> >, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens...>::template tuple_element_type<C> >::type;
                typedef typename Prepend<std::vector<TriggerDatum<typename TokenType::type, AllTokens...> >, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::tuple_type >::type tuple_type;
                template <char C> using index = typename std::conditional<C == Tag, Index<I>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens...>::template index<C> >::type;
            };
            
            template <char Tag, typename ... Tokens, unsigned int I, typename ... AllTokens>
            struct TriggerData<std::tuple<Trigger<Tag, VoidToken>, Tokens...>, I, AllTokens...>
            {
                template <char C> using token_type = typename std::conditional<C == Tag, void, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template token_type<C> >::type;
                template <char C> using datum_type = typename std::conditional<C == Tag, TriggerDatum<void, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template datum_type<C> >::type;
                template <char C> using tuple_element_type = typename std::conditional<C == Tag, TriggerDatum<void, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template tuple_element_type<C> >::type;
                typedef typename Prepend<TriggerDatum<void, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::tuple_type >::type tuple_type;
                template <char C> using index = typename std::conditional<C == Tag, Index<I>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens...>::template index<C> >::type;
            };
      
            template <typename TokenType, typename ... Tokens, unsigned int I, typename ... AllTokens>
            struct TriggerData<std::tuple<VertexLabel<TokenType>, Tokens...>, I, AllTokens... >
            {
                template <char C> using token_type = typename std::conditional<C == ' ', typename TokenType::type, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template token_type<C> >::type;
                template <char C> using datum_type = typename std::conditional<C == ' ', TriggerDatum<typename TokenType::type, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template datum_type<C> >::type;
                template <char C> using tuple_element_type = typename std::conditional<C == ' ', TriggerDatum<typename TokenType::type, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template tuple_element_type<C> >::type;
                typedef typename Prepend<TriggerDatum<typename TokenType::type, AllTokens...>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::tuple_type >::type tuple_type;
                template <char C> using index = typename std::conditional<C == ' ', Index<I>, typename TriggerData<std::tuple<Tokens...>, I+1, AllTokens... >::template index<C> >::type;
            };
            
            template <typename T, typename ... Tokens, unsigned int I, typename ... AllTokens>
            struct TriggerData<std::tuple<T, Tokens...>, I, AllTokens...>
            {
                template <char C> using token_type = typename TriggerData<std::tuple<Tokens...>, I, AllTokens... >::template token_type<C>;
                template <char C> using datum_type = typename TriggerData<std::tuple<Tokens...>, I, AllTokens... >::template datum_type<C>;
                template <char C> using tuple_element_type = typename TriggerData<std::tuple<Tokens...>, I, AllTokens... >::template tuple_element_type<C>;
                typedef typename TriggerData<std::tuple<Tokens...>, I, AllTokens... >::tuple_type tuple_type;
                template <char C> using index = typename TriggerData<std::tuple<Tokens...>, I, AllTokens...>::template index<C>;
            };
            
            //////////////////////
            
            template <typename T, typename ... Tokens>
            struct TriggerDatum
            {
                TriggerDatum()
                {
                }
                
                template <char OpeningTag>
                typename Detail::TriggerCommentData<std::tuple<Tokens...> >::template tuple_element_type<OpeningTag> & comments()
                {
                    return std::get<Detail::TriggerCommentData<std::tuple<Tokens...> >::template index<OpeningTag>::value>(commentData);
                }
                
                template <char OpeningTag>
                void addComment(typename Detail::TriggerCommentData<std::tuple<Tokens...> >::template token_type<OpeningTag> && comment)
                {
                    std::get<Detail::TriggerCommentData<std::tuple<Tokens...> >::template index<OpeningTag>::value>(commentData).emplace_back(std::move(comment));
                }
                
                const T & value() const
                {
                    return val;
                }
                
                void setValue(T && value) 
                {
                    val=std::move(value);
                }  
                
                unsigned int pos() const
                {
                    return p;
                }
                
                void setPos(unsigned int i)
                {
                    p = i;
                }
                
            private:
                
                T val;
                unsigned int p;
                typename TriggerCommentData<std::tuple<Tokens...> >::tuple_type commentData;
            };
            
            template <typename ... Tokens>
            struct TriggerDatum<void, Tokens...>
            {
                TriggerDatum()
                {
                }
                
                template <char OpeningTag>
                typename Detail::TriggerCommentData<std::tuple<Tokens...> >::template tuple_element_type<OpeningTag> & comments()
                {
                    return std::get<Detail::TriggerCommentData<std::tuple<Tokens...> >::template index<OpeningTag>::value>(commentData);
                }
                
                template <char OpeningTag>
                void addComment(typename Detail::TriggerCommentData<std::tuple<Tokens...> >::template token_type<OpeningTag> && comment)
                {
                    std::get<Detail::TriggerCommentData<std::tuple<Tokens...> >::template index<OpeningTag>::value>(commentData).emplace_back(std::move(comment));
                }
                
                unsigned int pos() const
                {
                    return p;
                }
                
                void setPos(unsigned int i)
                {
                    p = i;
                }
                
            private:
                
                unsigned int p;
                typename TriggerCommentData<std::tuple<Tokens...> >::tuple_type commentData;
            };
            
            
            template <typename ... Tokens>
            struct ElementData
            {                
                template <char Tag>
                typename Detail::TriggerData<std::tuple<Tokens...>, 0, Tokens...>::template tuple_element_type<Tag> & triggerData()
                {
                    return std::get<Detail::TriggerData<std::tuple<Tokens...>, 0, Tokens...>::template index<Tag>::value>(data);
                }
            
                template <char Tag>
                typename Detail::TriggerData<std::tuple<Tokens...>, 0, Tokens...>::template datum_type<Tag> & activeTriggerData()
                {
                    return getActiveTriggerData(std::get<Detail::TriggerData<std::tuple<Tokens...>, 0, Tokens...>::template index<Tag>::value>(data));
                }
                
                template <char Tag>
                void addTrigger(unsigned int pos)
                {
                    addTriggerData<Tag> (std::get<Detail::TriggerData<std::tuple<Tokens...>, 0, Tokens...>::template index<Tag>::value>(data));
                    activeTriggerData<Tag>().setPos(pos);
                }

            private:
                
                template <typename T>
                T & getActiveTriggerData(T & t)
                {
                    return t;
                }
                
                template <typename T, typename Alloc>
                T & getActiveTriggerData(std::vector<T, Alloc> & t)
                {
                    return t.back();
                }
                
                template <char Tag, typename T>
                void addTriggerData(T & t) 
                {
                }
                
                template <char Tag, typename T, typename Alloc>
                void addTriggerData(std::vector<T, Alloc> & t) 
                {
                    t.emplace_back(std::move(typename TriggerData<std::tuple<Tokens...>, 0, Tokens...>::template datum_type<Tag>()));
                }
                
                typename TriggerData<std::tuple<Tokens...>, 0, Tokens...>::tuple_type data;
            };
        }
    }
}

#endif // GULO_NEWICK_ELEMENTDATA_HPP