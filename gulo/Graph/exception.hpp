#ifndef GULO_GRAPH_EXCEPTION_HPP
#define GULO_GRAPH_EXCEPTION_HPP

#include <exception>
#include <string>

namespace Gulo
{
    class BadGraphElementException: public std::exception
    {
    public:

        BadGraphElementException(const std::string & message)
        : msg(message)
        {
        }

        virtual const char* what() const throw()
        {
            std::string result = "Bad graph element";
            if (msg.size() > 0) {
                result += ": " + msg;
            }
            return result.c_str();
        }

    private:

        std::string msg;
    };

    class BadGraphIsomorphismException: public std::exception
    {
    public:

        BadGraphIsomorphismException(const std::string & message)
        : msg(message)
        {
        }

        virtual const char* what() const throw()
        {
            std::string result = "Bad graph isomorphism";
            if (msg.size() > 0) {
                result += ": " + msg;
            }
            return result.c_str();
        }

    private:

        std::string msg;
    };
}

#endif // GULO_GRAPH_EXCEPTION_HPP