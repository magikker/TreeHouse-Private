#ifndef GULO_NEWICK_COMMON_HPP
#define GULO_NEWICK_COMMON_HPP

#include <climits>
#include <array>

namespace Gulo
{
    namespace Newick
    {
        template <char C, typename TokenType>
        class Trigger;
        
        template <char OpeningTag, char ClosingTag, typename TokenType>
        class Comment;
        
        template <typename TokenType>
        struct VertexLabel;
        
        template <template <typename...> class ProcessorType>
        struct Processor;
        
        struct ParseResult
        {
            ParseResult(const std::string & newickString) 
            : hasErr(false),
                errPos(0),
                newick(newickString)
            {
            }
            
            void setError(const std::string & message, unsigned int pos)
            {
                hasErr = true;
                errMessage = message;
                errPos = pos;
                static const int len = 60;
                std::string n = newick;
                n = n.substr((pos > len/2) ? pos - len/2 : 0, len);
                if (pos > len/2) 
                    pos = len/2;
                if (n.size() > 0) {
                    int pregap = (pos > n.size()) ? n.size() : pos;
                    errSnippet = n + "\n" + std::string(pregap, ' ') + "^";
                }
            }
            
            bool hasError() const
            {
                return hasErr;
            }
            
            const std::string & errorMessage() const
            {
                return errMessage;
            }
            
            const std::string & errorSnippet() const
            {
                return errSnippet;
            }
            
            unsigned int errorPos() const
            {
                return errPos;
            }
            
        private:
            
            bool hasErr;
            std::string errMessage;
            unsigned int errPos;
            std::string newick;
            std::string errSnippet;
        };
        
        struct CharacterValidator
        {
            void reserve(char ch) 
            {
                arr[ch-CHAR_MIN] = 1;
            }
            
            bool isReserved(char ch) const
            {
                return arr[ch-CHAR_MIN];
            }
            
            std::array<char, CHAR_MAX-CHAR_MIN+1> arr{};
        };
        
        struct VoidToken
        {
        };
        
        namespace Detail
        {
            template <typename T, typename VariadicContainer>
            struct Prepend;
            
            template <template <typename ...> class VariadicContainerType, typename T, typename ... T2>
            struct Prepend<T, VariadicContainerType<T2...> >
            {
                typedef VariadicContainerType<T, T2...> type;
            };
            
            template <unsigned int I>
            struct Index
            {
                static const unsigned int value = I;
            };
        }
    }
}

#endif // GULO_NEWICK_COMMON_HPP