using Nemo
Q,x = PolynomialRing(QQ,"x")
q=2#the power of the geometric series
function poly_ceil(f,h)
	r=zero(Q)
	for i=0:h-1
		r+=coeff(f,i)*x^i
	end
	return r
end
function poly_floor(f,l)
	r=zero(Q)
	for i=0:degree(f)-l
		r+=coeff(f,i+l)*x^i
	end
	return r
end
function t_kar_mul(n,a,c)
	m = degree(a)
	if(degree(c)> m + n)
		return "error"
	end
	if n==0
		r=zero(Q)
		for k=0:m
			r+=coeff(a,k)*coeff(c,k)
		end
		return r
	end
	if m<=0
		return coeff(a,0)*c*one(Q)
	end
	mu = div(m,2)+1
	vu = div(n,2)+1
	h = max(mu,vu)
	rho = degree(poly_ceil(a,h))
	if mu > n
		u = t_kar_mul(n,poly_ceil(a,mu),poly_ceil(c,n+rho))
		v = t_kar_mul(n,poly_floor(a,mu),poly_floor(c,mu))
		return u+v * one(Q)
	end
	if vu > m
		u = t_kar_mul(vu-1,a,poly_ceil(c,m+vu))
		v = t_kar_mul(n-vu,a,poly_floor(c,vu))
		return u + v * one(Q)
	end
	lambda = degree(poly_ceil(a,h) + poly_floor(a,h))
	r = t_kar_mul(h-1,poly_ceil(a,h),poly_ceil(c,rho+h)-poly_floor(poly_ceil(c,rho+2*h),h))
	s = t_kar_mul(n-h,poly_floor(a,h),poly_floor(c,2*h)-poly_floor(poly_ceil(c,m+n-h+1),h))
	t = t_kar_mul(h-1,poly_ceil(a,h) + poly_floor(a,h), poly_floor(poly_ceil(c,lambda + 2*h),h))
	return r + t + x^h*(s+poly_ceil(t,n-h+1))
end
#only works for n = 2^k, k \in Z right now i think
function fast_mat_mult(A,B)
	n = size(A)[1]
	if n==1
		return [A[1,1] * B[1,1]]
	end
	half_n = div(n,2)
	A11 = A[1:half_n,1:half_n]
	A12 = A[1:half_n,half_n+1:n]
	A21 = A[half_n+1:n,1:half_n]
	A22 = A[half_n+1:n,half_n+1:n]
	B11 = B[1:half_n,1:half_n]
	B12 = B[1:half_n,half_n+1:n]
	B21 = B[half_n+1:n,1:half_n]
	B22 = B[half_n+1:n,half_n+1:n]
	
	S1 = A21+A22
	S2 = S1-A11
	S3 = A11-A21
	S4 = A12-S2
	T1 = B12-B11
	T2 = B22-T1
	T3 = B22-B12
	T4 = T2-B21
	
	P1 = fast_mat_mult(A11,B11)
	P2 = fast_mat_mult(A12,B21)
	P3 = fast_mat_mult(S4,B22)
	P4 = fast_mat_mult(A22,T4)
	P5 = fast_mat_mult(S1,T1)
	P6 = fast_mat_mult(S2,T2)
	P7 = fast_mat_mult(S3,T3)

	U1 = P1 + P2
	U2 = P1 + P6
	U3 = U2 + P7
	U4 = U2 + P5
	U5 = U4 + P3
	U6 = U3 - P4
	U7 = U3 + P5
	
	return [U1 U5; U6 U7]
end
function newton_to_monomial_geom(fs)
	n = length(fs)
	qs=Array{Nemo.fmpq}(n)
	vs=Array{Nemo.fmpq}(n)
	us=Array{Nemo.fmpq}(n)
	ws=Array{Nemo.fmpq}(n)
	zs=Array{Nemo.fmpq}(n)
	qs[1]=1
	vs[1]=fs[1]
	us[1]=1
	ws[1]=fs[1]
	zs[1]=1
	qi=coeff(one(Q),0)
	negi=-1
	for i=2:n
		qs[i]=qs[i-1]*qi
		qi = q*qi
		vs[i]=(negi)*fs[i]*qs[i]
		us[i]=div(us[i-1]*qi,(1-qi))
		ws[i]=div(vs[i],us[i])
		zs[i]=negi*div(us[i],qs[i])
		negi*=-1
	end
	G = t_kar_mul(n,reduce(+,0,[x^i*us[i+1] for i in 0:n-1]),reduce(+,0,[x^i*ws[i+1] for i in 0:n-1]))
	return [zs[i+1]*coeff(G,i) for i in 0:n-1]
end
function newton_interp_geom(vs)
	n = length(vs)
	qs = Array{Nemo.fmpq}(n)
	us = Array{Nemo.fmpq}(n)
	ws = Array{Nemo.fmpq}(n)
	qs[1]=1
	us[1]=1
	ws[1]=vs[1]
	qi = 1
	for i=2:n
		qs[i] = qs[i-1]*qi
		qi*=q
		us[i] = us[i-1]*(qi-1)
		ws[i] = div(vs[i],us[i])
	end
	f = zero(Q)
	g = zero(Q)
	neg = 1
	for i=0:n-1
		f += x^i * ws[i+1]
		g += neg * x^i * div(qs[i+1],us[i+1])
		neg *=-1
	end
	G =	f*g 
	return [div(coeff(G,i),qs[i+1]) for i in 0:n-1]
end

#takes an array of a partioned polynomial and transforms it into a matrix
function fs_to_mat(fs,n,m)
	A = Array{typeof(coeff(fs[1],1))}(m,n)
	for i=1:m
		for j=1:n
			A[i,j]=coeff(fs[i],j-1)
		end
	end
	return A
end

function fast_comp_mod(f,g,h)
	n = degree(f)
	N=2*n
	m=convert(Int64,ceil(n^(1/2)))
	G = typeof(g)[]
	for i=1:m
		gi = g % x^(m) 
		g = divexact((g- gi),x^(i*m))
		push!(G, gi)
	end
	Λ = Array{Nemo.fmpq}(m,N)
	H = Array{Nemo.fmpq_poly}(m)
	H[1] = one(Q)
	for i=2:m
		H[i]=(H[i-1]*h) % f
	end
	#A = fs_to_mat(H,n,m)
	A = Array{Nemo.fmpq}(m,N)
	#qs are powers of q
	for i=1:m
		qj=1
		for j=1:N
			A[i,j]=H[i](qj)
			qj = q*qj
		end
	end
	B = fs_to_mat(G,m,m)
#	BA = fast_mat_mult(B,A)
	BA = B*A
	b=0
	hm = H[m]*h
	Λi = one(Q)
	for i=1:m
		qi=1
		for j=1:N
			Λ[i,j]=Λi(qi)
			qi*=q
		end
		Λi=(Λi*hm)%f
	end
	r = reduce(+, 0, [BA[i,:] .* Λ[i,:] for i in 1:m])
	coeffs =  newton_to_monomial_geom(newton_interp_geom(r))
	rho = zero(Q)
	for i=0:length(coeffs)-1
		rho += x^i * coeffs[i+1]
	end
	return rho%f
end

f = 3+2x^2+x^4+x^10
g = 1+x+x^3
h = 1+x^2+x^4
@time println(fast_comp_mod(f,g,h))
@time println(g(h)%f)
