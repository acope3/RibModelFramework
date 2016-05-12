#ifndef Utility_H
#define Utility_H
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


// This template handles general printing between C++ and R
// Returns 0 if no errors in formatting detected.
inline int my_print(const char *s)
{
    int rv = 0; // By default, assume success

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
               //throw std::runtime_error("invalid format string: missing arguments");
               rv = 1;
        }
#ifndef STANDALONE
        Rcpp::Rcout << *s++;
#else
        std::cout << *s++;
#endif
    }

    return rv;
}

template<typename T, typename... Args>
inline int my_print(const char *s, T value, Args... args)
{
    int rv = 0;

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
            {
#ifndef STANDALONE
                Rcpp::Rcout << value;
#else
                std::cout << value;
#endif
                rv = my_print(s + 1, args...); // call even when *s == 0 to detect extra arguments
                return rv;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcout << *s++;
#else
        std::cout << *s++;
#endif
    }
    //throw std::logic_error("extra arguments provided to my_print");
    return 1;
}

// This template handles error printing between C++ and R
// Returns 0 if no errors in formatting detected.
inline int my_printError(const char *s)
{
    int rv = 0;

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
                //throw std::runtime_error("invalid format string: missing arguments");
                rv = 1;
        }
#ifndef STANDALONE
        Rcpp::Rcerr << *s++;
#else
        std::cerr << *s++;
#endif
    }

    return rv;
}

template<typename T, typename... Args>
inline int my_printError(const char *s, T value, Args... args)
{
    int rv = 0;

    while (*s)
    {
        if (*s == '%')
        {
            if (*(s + 1) == '%')
                ++s;
            else
            {
#ifndef STANDALONE
                Rcpp::Rcerr << value;
#else
                std::cerr << value;
#endif
                rv = my_printError(s + 1, args...); // call even when *s == 0 to detect extra arguments
                return rv;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcerr << *s++;
#else
        std::cerr << *s++;
#endif
    }
    //throw std::logic_error("extra arguments provided to my_print");
    return 1;
}

//Blank header
#endif // Utility_H
