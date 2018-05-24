using Nemo
#Z,x = PolynomialRing(ZZ,"x")
Zp = ResidueRing(ZZ, 109)
Z, x = PolynomialRing(Zp, "x")
function poly_ceil(f,h)
	return truncate(f,h)
end
function poly_floor(f,l)
	return shift_right(f,l)
end
function t_kar_mul(n,a,c)
	m = degree(a)
	#println("n is $n m is $m deg(c) is $(degree(c))")
	#println("c is $c")
	#println("a is $a")
	if(degree(c)> m + n)
		return "error"
	end
	if n==0
		r=zero(Z)
		for k=0:m
			r+=coeff(a,k)*coeff(c,k)
		end
		return r
	end
	if m<=0
		return coeff(a,0)*c*one(Z)
	end
	mu = div(m,2)+1
	vu = div(n,2)+1
	h = max(mu,vu)
	rho = h-1
	#rho = h
	if mu > n
		u = t_kar_mul(n,poly_ceil(a,mu),poly_ceil(c,n+rho))
		v = t_kar_mul(n,poly_floor(a,mu),poly_floor(c,mu))
		return u+v * one(Z)
	end
	if vu > m
		u = t_kar_mul(vu-1,a,poly_ceil(c,m+vu))
		v = t_kar_mul(n-vu,a,poly_floor(c,vu))
		return u + v * one(Z)
	end
	lambda = degree(poly_ceil(a,h) + poly_floor(a,h))
	#println("rho $rho h $h")
	#println("poly $(poly_ceil(c,rho+h))")
	#nprintln("poly2 $(poly_floor(poly_ceil(c,rho+2*h),h))")
	#println("r rec")
	r = t_kar_mul(h-1,poly_ceil(a,h),poly_ceil(c,rho+h)-poly_floor(poly_ceil(c,rho+2*h),h))
	#println("s rec");
	s = t_kar_mul(n-h,poly_floor(a,h),poly_floor(c,2*h)-poly_floor(poly_ceil(c,m+n-h+1),h))
	#println("rec 2")
	#println("t rec");
	t = t_kar_mul(h-1,poly_ceil(a,h) + poly_floor(a,h), poly_floor(poly_ceil(c,lambda + 2*h),h))
	#println("rec 3")
	println("p4 $(s+poly_ceil(t,n-h+1))")
	println("p6 $(x^h*(s+poly_ceil(t,n-h+1)))")
	println("p5 $(t+x^h*(s+poly_ceil(t,n-h+1)))")
	return r + t+x^h*(s+poly_ceil(t,n-h+1))
end
function newton_to_monomial_geom(fs,q)
	n = length(fs)
	qs=Array{Nemo.nmod}(n)
	vs=Array{Nemo.nmod}(n)
	us=Array{Nemo.nmod}(n)
	ws=Array{Nemo.nmod}(n)
	zs=Array{Nemo.nmod}(n)
	qs[1]=lead(one(Z))
	vs[1]=fs[1]
	us[1]=lead(one(Z))
	ws[1]=fs[1]
	zs[1]=lead(one(Z))
	qi=coeff(one(Z),0)
	negi=-1
	for i=2:n
		qs[i]=qs[i-1]*qi
		qi = q*qi
		vs[i]=(negi)*fs[i]*qs[i]
		us[i]=us[i-1]*qi*inv(1-qi)
		ws[i]=vs[i]*inv(us[i])
		zs[i]=negi*us[i]*inv(qs[i])
		negi*=-1
	end
	#println(reduce(+,0,[x^i*us[i+1]for i in 0:n-1]))
	#println(reduce(+,0,[x^i*ws[i+1]for i in 0:n-1]))
	G = t_kar_mul(n,reduce(+,0,[x^i*us[i+1] for i in 0:n-1]),reduce(+,0,[x^i*ws[i+1] for i in 0:n-1]))
	println("G is $G")
	println("n is $n")
	println("q is $qs")
	r = [zs[i+1]*coeff(G,i) for i in 0:n-1]
	println(r)
	return r
end
function newton_interp_geom(vs,q)
	n = length(vs)
	qs = Array{Any}(n)
	us = Array{Any}(n)
	ws = Array{Any}(n)
	qs[1]=lead(one(Z))
	us[1]=lead(one(Z))
	ws[1]=vs[1]
	qi = 1
	for i=2:n
		qs[i] = qs[i-1]*qi
		qi*=q
		us[i] = us[i-1]*(qi-1)
		ws[i] = vs[i]*inv(us[i])
	end
	f = zero(Z)
	g = zero(Z)
	neg = 1
	for i=0:n-1
		f += x^i * ws[i+1]
		g += neg * x^i * qs[i+1]*inv(us[i+1])
		neg *=-1
	end
	G =	f*g 
	ret= [coeff(G,i)*inv(qs[i+1]) for i in 0:n-1]
	return ret
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
	m=convert(Int64,ceil(n^(1/2)))
	N=2*n
	p = convert(Int64,modulus(Z))
	q=2
	for qi = 2:p
		k = 1
		qk = qi
		while (qk != 1 && k < N)
			qk=(qk*qi)%p
			k+=1
		end
		if k >= N
			q = qi
			break
		end
	end
	println(q)
	G = typeof(g)[]
	for i=1:m
		gi = g % x^(m) 
		g = divexact((g- gi),x^(i*m))
		push!(G, gi)
	end
	Λ = Array{Any}(m,N)
	H = Array{Any}(m)
	H[1] = one(Z)
	for i=2:m
		H[i]=(H[i-1]*h) % f
	end
	#A = fs_to_mat(H,n,m)
	A = Array{Any}(m,N)
	#qs are powers of q
	for i=1:m
		qj=1
		for j=1:N
			A[i,j]=H[i](qj)
			qj = q*qj
		end
	end
	B = fs_to_mat(G,m,m)
	BA = B*A
	b=0
	hm = H[m]*h
	Λi = one(Z)
	for i=1:m
		qi=lead(one(Z))
		for j=1:N
			Λ[i,j]=Λi(qi)
			qi*=q
		end
		Λi=(Λi*hm)%f
	end
	r = reduce(+, 0, [BA[i,:] .* Λ[i,:] for i in 1:m])
	coeffs =  newton_to_monomial_geom(newton_interp_geom(r,q),q)
	rho = zero(Z)
	for i=0:length(coeffs)-1
		rho += x^i * coeffs[i+1]
	end
	return rho%f
end

f = 3+2x^2+x^4
g = 1+x+x^3
h = 2+x^2
@time println(fast_comp_mod(f,g,h))
@time println(g(h)%f)
