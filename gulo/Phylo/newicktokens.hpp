#ifndef GULO_NEWICK_COMMENT_HPP
#define GULO_NEWICK_COMMENT_HPP

#include <string>
#include <cstdlib>
#include <cctype>
#include <climits>
#include "newickcommon.hpp"

namespace Gulo
{
    namespace Newick
    {
        struct NestedComment
        {      
            typedef std::string type;
                      
            static type process(const std::string & str, unsigned int & i, char openingTag, char closingTag, ParseResult & result, const CharacterValidator &)
            {
                if (str[i] != openingTag) {
                    result.setError("Expected " + std::string(1, openingTag) + " at the beginning of a comment", i);
                    return "";
                }
                int commentCount = 1;
                unsigned int start = i;
                ++i;
                while (i < str.size() && commentCount > 0) {
                    if (str[i] == openingTag) {
                        ++commentCount;
                    } else if (str[i] == closingTag) {
                        --commentCount;
                    }
                    ++i;
                }
                if (i == str.size() || commentCount > 0) {
                    result.setError("Comment was opened but not closed", start);
                } else if (i - start - 2 != 0) {
                    return str.substr(start + 1, i - start - 2);
                }
                return "";
            }
        };

        struct UnnestedComment
        {      
            typedef std::string type;
            
            static type process(const std::string & str, unsigned int & i, char openingTag, char closingTag, ParseResult & result, const CharacterValidator &)
            {
                
                if (str[i] != openingTag) {
                    result.setError("Expected " + std::string(1, openingTag) + " at the beginning of a comment", i);
                    return "";
                }
                unsigned int start = i;
                ++i;
                while (i < str.size() && str[i] != closingTag && str[i] != openingTag) {
                    ++i;
                }
                if (i == str.size()) {
                    result.setError("Comment was opened but not closed", start);
                    return "";
                } else if (str[i] == openingTag) {
                    result.setError("Nested comments are not permitted", i);
                    return "";
                } else if (i - start - 1 != 0) {
                    ++i;
                    return str.substr(start + 1, i - start - 2);
                }
                return "";
            } 
        };
        
        struct UnquotedToken
        {  
            typedef std::string type;
            
            static type process(const std::string & str, unsigned int & i, ParseResult & result, const CharacterValidator & validator)
            {
                
                char ch = str[i];
                while (i < str.size() && std::isspace(ch)) {
                    ch = str[++i];
                }
                unsigned int start = i;
                while (i < str.size() && !validator.isReserved(ch) && !std::isspace(ch) && str[i] != '\'') {
                    ch = str[++i];
                }          
                if (i < str.size() && ch == '\'') {
                    result.setError("Unexpected quotation mark in unquoted token", i);
                    return "";
                }
                if (i-start > 0) {
                    std::string label = str.substr(start, i-start);
                    for (auto & c : label) {
                        if (c == '_') {
                            c = ' ';
                        }
                    }
                    
                    return label;
                }
                return "";
            }
        };
        
        
        struct QuotedToken
        {  
            typedef std::string type;
            
            static type process(const std::string & str, unsigned int & i, ParseResult & result, const CharacterValidator & validator)
            {             
                
                char ch = str[i];
                while (i < str.size() && std::isspace(ch)) {
                    ch = str[++i];
                }
                unsigned int start = i;
                if (i == str.size() || ch != '\'') {
                    result.setError("Expected quoted token to begin with an apostrophe", i);
                    return "";
                }
                ++i;
                while (i < str.size()-1 && str[i] == '\'' && str[i+1] == '\'') {
                    i += 2;
                }
                while ((i < str.size()) && str[i] != '\'') {
                    ++i;
                    while (i < str.size()-1 && str[i] == '\'' && str[i+1] == '\'') {
                        i += 2;
                    }
                }
                if (i == str.size()) {
                    result.setError("Quoted token was opened but not closed", start);
                    return "";
                }            
                ++i;
                if (i > start + 2) {
                    std::string label = str.substr(start+1, i-start-2);
                    
                    return label;
                }
                return "";
            }
        };
        
        struct UnquotedTokenWithApostrophes
        {  
            typedef std::string type;
            
            static type process(const std::string & str, unsigned int & i, ParseResult & result, const CharacterValidator & validator)
            {
                char ch = str[i];
                while (i < str.size() && std::isspace(ch)) {
                    ch = str[++i];
                }
                unsigned int start = i;
                while (i < str.size() && !validator.isReserved(ch)) {
                    ch = str[++i];
                }  
                if (i != start) {
                    std::string label = str.substr(start, i-start);
                    for (auto & c : label) {
                        if (c == '_') {
                            c = ' ';
                        }
                    }
                    return label;
                }
                return "";
            }
        };
        
        struct QuotedOrUnquotedToken
        {  
            typedef std::string type;
            
            static type process(const std::string & str, unsigned int & i, ParseResult & result, const CharacterValidator & validator)
            {             
                while (i < str.size() && std::isspace(str[i])) {
                    ++i;
                }
                if (i < str.size() && str[i] == '\'') { 
                    return QuotedToken::process(str, i, result, validator);
                } else {
                    return UnquotedToken::process(str, i, result, validator);
                }
    
            }
        };
        
        struct IntegerToken
        {  
            typedef int type;
            
            static type process(const std::string & str, unsigned int & i, ParseResult & result, const CharacterValidator & validator)
            {
                
                while (i < str.size() && std::isspace(str[i])) {
                    ++i;
                }
                unsigned int start = i;
                if (i < str.size() && (str[i] == '-' || str[i] == '+')) {
                    ++i;
                }
                while (i < str.size() && std::isdigit(str[i])) {
                    ++i;
                }   
                if (i < str.size() && !validator.isReserved(str[i]) && !std::isspace(str[i])) {
                    result.setError("Unexpected character in integer token", i);
                } else if (i > start) {
                    std::string token = str.substr(start, i-start);
                    char * e = 0;
                    long int d = std::strtol(token.c_str(), &e, 10);
                    if (!(e != 0 && *e == 0)) {
                        result.setError("Could not convert token to integer type", start);  
                    } else if (errno == ERANGE && d > INT_MAX) {
                        result.setError("Integer overflow", start);
                    } else if (errno == ERANGE && d < INT_MIN) {
                        result.setError("Integer underflow", start);                
                    } else {
                        
                        return (int)d;
                    }
                }
                return 0;
            }
        };
        
        
        struct RealToken
        {  
            typedef double type;
            
            static type process(const std::string & str, unsigned int & i, ParseResult & result, const CharacterValidator & validator)
            {
                
               char ch = str[i];
                while (i < str.size() && std::isspace(ch)) {
                    ch = str[++i];
                }
                unsigned int start = i;
                if (i < str.size() && (ch == '-' || ch == '+')) {
                    ch = str[++i];
                }
                while (i < str.size() && std::isdigit(ch)) {
                    ch = str[++i];
                }   
                if (i < str.size() && ch == '.') {
                    ch = str[++i];
                    while (i < str.size() && std::isdigit(ch)) {
                        ch = str[++i];
                    } 
                }
                if (i < str.size() && (ch == 'e' || ch == 'E')) {
                    ch = str[++i];
                    if (i < str.size() && (ch == '+' || ch == '-')) {
                        ch = str[++i];
                    }
                    if (i == str.size() || !std::isdigit(ch)) {
                        result.setError("Expected digit after exponent in real token", i);
                        return 0;
                    }
                    ch = str[++i];
                    while (i < str.size() && std::isdigit(ch)) {
                        ch = str[++i];
                    }
                }            
                if (i < str.size() && !validator.isReserved(ch) && !std::isspace(ch)) {
                    result.setError("Unexpected character in real token", i);
                } else if (i != start) {
                    std::string token = str.substr(start, i-start);
                    char * e = 0;
                    double d = std::strtod(token.c_str(), &e);
                    if (e == 0 || *e != 0) {
                        result.setError("Could not convert token to real type", start);
                    } else {
                        
                        return d;
                    }
                }
                return 0;
            }
        };
    }
}
    
#endif // GULO_NEWICK_COMMENT_HPP