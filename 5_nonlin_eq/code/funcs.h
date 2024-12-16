#pragma once

#include<iostream>
#include <array>
#include <type_traits>
#include <tuple>

struct IterationException {};

/**
    Ðåøàåò óðàâíåíèå Êåïëåðà ìåòîäîì Íüþòîíà
    * ecc - ýêñöåíòðèñèòåò, ïðèíàäëåæèò (0, 1)
    * meanAnomaly - ñðåäíÿÿ àíîìàëèÿ, Ì (ðàäèàíû)
    * maxIter - ìàêñèìàëüíîå êîëè÷åñòâî èòåðàöèé
    * tol - òî÷íîñòü, ñ êîòîðîé íóæíî îòûñêàòü ðåøåíèå

    Ðåêîìåíäóåìîå ïîâåäåíèå. Åñëè ðåøåíèå íå íàøëîñü çà maxIter èòåðàöèé - âûáðàñûâàòü èñêëþ÷åíèå.
    Åñëè ïðèáëèæåíèÿ ê ðåøåíèþ ìåæäó èòåðàöèÿìè ìåíÿþòñÿ íå áîëåå, ÷åì íà tol, òî ðåøåíèå äîñòèãíóòî.
**/
double keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol)
{
    double E = meanAnomaly + ecc * std::sin(meanAnomaly) /
        (1 - ecc * std::cos(meanAnomaly));
    double delta = std::abs(meanAnomaly - E);

    unsigned int i;
    for (i = 1; i < maxIter and delta > tol; i++) {
        delta = E;
        E = E - (E - ecc * std::sin(E) - meanAnomaly) /
            (1 - ecc * std::cos(E));
        delta = std::abs(delta - E);
    }
    
    //if (i == maxIter) { throw IterationException(); };

    return E;
};



template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(
    const Callable& func,                                             
    const RealType& tau,                                              
    const typename ArgumentGetter<Callable>::Argument& initialGuess,  
    const unsigned int nIteration                                     
) 
{
    double x = initialGuess;

    for (int i = 0; i < nIteration; i++) {
        x += tau * func(x);
    }

    return x;

};
