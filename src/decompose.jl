###
### Series as linear functionals on polynomials. 
###
### Decomposition of series from their first moments.
### Compute the corresponding weighted sum of Dirac measures.
###
### Bernard Mourrain
###


export decompose

#------------------------------------------------------------------------
function ratio(V1,V0)
   i = 1;
   while V0[i] == 0 
      i+=1;
   end;
   V1[i]/V0[i];
end

#------------------------------------------------------------------------
function numrank(S, eps)
  i :: Int = 1;
  while i<= length(S) && S[i]/S[1] > eps
    i+= 1;
  end
  i-1;
end

#------------------------------------------------------------------------
function decompose(sigma :: Series{C,Mon}, U, S, V, B0, B1, X, r::Int64) where {C,Mon}

    n = length(X)
    Un0 = transpose(U[1,1:r])
    Un1 = V[1,1:r]

    Sinv = diagm([one(C)/S[i] for i in 1:r])
    
    M = []
    for i in 1:n
        H = hankel(sigma, B0, B1*X[i])
    	push!(M, Sinv*(ctranspose(U[:,1:r])*H*V[:,1:r]))
    end
    
    M0 = M[1]
    for i in 2:n
       M0 += M[i]*rand(Float64)
    end

    E = eigvecs(M0)

    #println("Eigenval=", eigvals(M0))

    Xi = fill(zero(E[1,1]),r,n)

    w = E \ Un1

    for i in 1:r
    	for j in 1:n
	    Xi[i,j] = ratio(M[j]*E[:,i], E[:,i])
	end
    end

    D = Un0*diagm([S[i] for i in 1:r])

    for i in 1:r
    	w[i] *=  (D*E[:,i])[1]
    end

    return w, Xi
end
#------------------------------------------------------------------------
"""
```
decompose(σ :: Series{C,M}, eps :: Float64 = 1.e-6)
```
Decompose the series ``σ`` as a weighted sum of exponentials.
Return ``ω``, ``Ξ`` where 
 - ``ω`` is the vector of weights,
 - ``Ξ`` is the matrix of frequency points, stored per row.
The list of monomials of degree ``\\leq {d-1 \\over 2}`` are used to construct 
the Hankel matrix, where ``d`` is the maximal degree of the moments in ``σ``.

The optional argument ``eps`` is the precision used to compute the numerical rank.
Its default value is ``1.e-6``.
"""
function decompose(sigma :: Series{C,M}, eps :: Float64 = 1.e-6,
                       B0 = monoms(variables(sigma), div(maxdegree(sigma)-1,2) ),
                       B1 = monoms(variables(sigma), div(maxdegree(sigma)-1,2) ),
                       X  = variables(sigma)) where {C,M}
    H = hankel(sigma, B0, B1)
    U, S, V = svd(H)       #H= U diag(S) V'
    r = numrank(S, eps)
    decompose(sigma, U,S,V, B0, B1, X, r)
end

#------------------------------------------------------------------------
"""
```
decompose(σ :: Series{C,M}, r::Int64)
```
Decompose the series ``σ`` assuming its rank is r.
"""
function decompose(sigma :: Series{C,M}, r :: Int64,
                       B0 = monoms(variables(sigma), div(maxdegree(sigma)-1,2) ),
                       B1 = monoms(variables(sigma), div(maxdegree(sigma)-1,2) ),
                       X  = variables(sigma)) where {C,M}
    H = hankel(sigma, B0, B1)
    U, S, V = svd(H)       #H= U diag(S) V'
    decompose(sigma, U,S,V, B0, B1, X, r)
end

