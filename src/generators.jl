###
### Series as linear functionals on polynomials. 
###
### Generating series from moment functions, which maps monomial
### exponents to values.
###
### Bernard Mourrain
###

export series, dual, moment, monoms

#-----------------------------------------------------------------------
"""
```
monoms(V, d::Int64) -> Vector{Monomial}
monoms(V, rg::UnitRangeInt64) -> Vector{Monomial}
```
List of all monomials in the variables V up to degree d of from degree d1 to d2, 
ordered by increasing degree.
"""
function monoms(V::Vector{Var}, rg::UnitRange{Int64}) where Var<:AbstractVariable
    l = [monomials(V,0)[1]]
    for i in 1:rg[end]
        append!(l,monomials(V,i))
    end
    l
end

function monoms(V::Vector{Var}, d ::Int64) where Var<:AbstractVariable
    if (d>0)
        monoms(V,0:d)
    else
        L = monoms(V, 0:-d)
        for i in 1:length(L)
            inv!(L[i])
        end
        L
    end
end

#-----------------------------------------------------------------------
"""
```
dual(p:Polynomial) -> Series{C,M}
```
Construct the series from the polynomial by replacing the monomial basis by its dual basis

Example
-------
```julia
@ring x1 x2
p = x1^2*x2+3*x2+2
s = dual(p)

# output 

dx1^2dx2 + 3dx2 + 2
```

"""
function dual(p::P) where P <: AbstractPolynomial
    C = coefficienttype(P)
    M = monomialtype(P)
    s = Series{C,M}()
    for t in terms(p)
        s[monomial(t)] = coefficient(t)
    end
    return s
end

#------------------------------------------------------------------------
""" 
```
series(f,L) -> Series{C,M}
```
Compute the generating series ``\\sum_{x^{α} \\in L} f(α) z^α``
for a function  ``f: \\mathbb{N}^n \\rightarrow T`` and a sequence L of monomials.
"""
function series(f::Function, L::AbstractVector)

    C = typeof(f(exponents(L[1])))
    M = eltype(L)
    res = series(f(exponents(L[1])), L[1])
   for m in L
         res[m] = f(exponents(m))
   end
   res
end;

#----------------------------------------------------------------------
function binom(d, alpha::Vector{Int64})
  r = binomial(d, alpha[1])
  for i in 2:length(alpha)
      d -= alpha[i-1]
      r *= binomial(d, alpha[i])
  end
  r
end

"""
```
series(p::Polynomial, d:: Int64) -> Series{C,M} 
```
Compute the series associated to the tensor p of degree d. 
C is the type of the coefficients of the polynomial p.
"""
function series(p::P, x = variables(p)[1], d:: Int64 = maxdegree(p)) where P<:AbstractPolynomial
    C = coefficienttype(P)
    M = monomialtype(P)
    s = Series{C,M}()
    pa = subs(p,x=>one(C))
    for t in pa
        s[monomial(t)] = coefficient(t)/binom(d,exponents(monomial(t)))
    end
    return s
end

#----------------------------------------------------------------------
"""
```
moment(w::Vector{C}, P::Matrix{C}) -> Vector{Int64} -> C
```
Compute the moment function ``α -> ∑_{i} ω_{i} P_{i}^α`` 
associated to the sequence P of r points of dimension n, which is a matrix 
of size r*n and the weights w.
"""
moment(w, P) = function(α::Vector{Int})
  res = 0
  for i in 1:size(P,1)
      m = 1
      for j in 1:length(α)
        m *= P[i,j]^α[j]
      end
      res+=m*w[i];
  end
  res
end

#----------------------------------------------------------------------
"""
```
moment(p::Polynomial, zeta::Vector{C}) -> Vector{Int64} -> C
```
Compute the moment function ``α \\rightarrow p(ζ^α)``.
"""
moment(p::P, zeta::Vector) where P<:AbstractPolynomial =  function(V::Vector{Int})
    U = [zeta[i]^V[i] for i in 1:length(V)]
    p(U)
end

"""
```
series(w:: Vector{C}, P::Matrix{C}, L::Vector{M}) -> Series{C,M}
```
Compute the series of the moment sequence ``∑_{i} ω_{i} P_{i}^α`` for ``α \\in L``.
"""
function series(w:: Vector{C}, P, L::Vector{M}) where {C,M}
   series(moment(w,P), L) 
end
    
#----------------------------------------------------------------------
"""
```
series(w:: Vector{C}, P::Matrix{C}, X, d::Int64) -> Series{C,M}
```
Compute the series of the moment sequence ``∑_i ω_{i} P_{i}^α`` for ``|α| \\leq d``.
"""
function series(w::Vector{C}, P::Matrix{C}, X, d::Int64) where C
    h = moment(w,P)
    L = monoms(X,d)
   series(h,L)
end


#----------------------------------------------------------------------
"""
```
series(p::Polynomial, zeta, X, d::Int64) -> Series{C,M}
```
Compute the series of moments ``p(ζ^α)`` for ``|α| \\leq d``.
"""
function series(p::P, zeta, X, d::Int64) where P<:AbstractPolynomial
    h = moment(p,zeta)
    L = monoms(X,d)
    series(h,L)
end
