
using QuadGK

function BlackScholes(S, T, K, r, σ)
    f(x) = exp(-0.5*x^2)
    N(d) = 1/sqrt(2*pi) * quadgk(f, -Inf, d)[1]
    d1 = ( log(S/K) + (r + 0.5 * σ^2)*(T))/(σ*sqrt(T))
    d2 = d1 - σ*sqrt(T)
    function call()
        return  S * N(d1) - K * exp(-r*(T))*N(d2)
    end
    function put()
        return K*exp(-r*(T))*N(-d2) - S*N(-d1)
    end
    return call(), put()
end

function BlackScholesDividend(S, T, K, r, σ, q)
    f(x) = exp(-0.5*x^2)
    N(d) = 1/sqrt(2*pi) * quadgk(f, -Inf, d)[1]
    d1 = ( log(S/K) + (r - q + 0.5 * σ^2)*(T))/(σ*sqrt(T))
    d2 = d1 - σ*sqrt(T)
    function call()
        return S*exp(-q*(T))*N(d1) - K*exp(-r*(T))*N(d2)    
    end
    function put()
        return K*exp(-r*(T))*N(-d2) - S*exp(-q*(T))*N(d1)
    end
    return call(), put()
end

BlackScholes(50,1,100,0.05,0.25)
