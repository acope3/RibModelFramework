#ifndef Utility_H
#define Utility_H
#include <iostream>


#ifndef STANDALONE
#include <Rcpp.h>
using namespace Rcpp;
#endif


// This template handles general printing between C++ and R
inline void my_print(const char *s)
{
    while (*s) {
        if (*s == '%') {
            if (*(s + 1) == '%') {
                ++s;
            }
            /*
               else {
               throw std::runtime_error("invalid format string: missing arguments");
               }
             */
        }
#ifndef STANDALONE
        Rcpp::Rcout << *s++;
#else
        std::cout << *s++;
#endif
    }
}

template<typename T, typename... Args>
inline void my_print(const char *s, T value, Args... args)
{
    while (*s) {
        if (*s == '%') {
            if (*(s + 1) == '%') {
                ++s;
            }
            else {
#ifndef STANDALONE
                Rcpp::Rcout << value;
#else
                std::cout << value;
#endif
                my_print(s + 1, args...); // call even when *s == 0 to detect extra arguments
                return;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcout << *s++;
#else
        std::cout << *s++;
#endif
    }
    //throw std::logic_error("extra arguments provided to printf");
}


// This template handles general printing between C++ and R
inline void my_printError(const char *s)
{
    while (*s) {
        if (*s == '%') {
            if (*(s + 1) == '%') {
                ++s;
            }
            /*
               else {
               throw std::runtime_error("invalid format string: missing arguments");
               }
             */
        }
#ifndef STANDALONE
        Rcpp::Rcerr << *s++;
#else
        std::cerr << *s++;
#endif
    }
}

template<typename T, typename... Args>
inline void my_printError(const char *s, T value, Args... args)
{
    while (*s) {
        if (*s == '%') {
            if (*(s + 1) == '%') {
                ++s;
            }
            else {
#ifndef STANDALONE
                Rcpp::Rcerr << value;
#else
                std::cerr << value;
#endif
                my_printError(s + 1, args...); // call even when *s == 0 to detect extra arguments
                return;
            }
        }
#ifndef STANDALONE
        Rcpp::Rcerr << *s++;
#else
        std::cerr << *s++;
#endif
    }
    //throw std::logic_error("extra arguments provided to printf");
}

//Blank header
#endif // Utility_H
