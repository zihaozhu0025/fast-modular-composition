using Nemo
Zp = ResidueRing(ZZ, 9001)
Z, x = PolynomialRing(Zp, "x")
function poly_ceil(f::Nemo.nmod_poly,h::Int64)
	return truncate(f,h)
end
function poly_floor(f::Nemo.nmod_poly,l::Int64)
	return shift_right(f,l)
end
function t_mul(n::Int64,a::Nemo.nmod_poly, c::Nemo.nmod_poly)
	m = degree(a)
	if( degree(c) > m+ n)
		println("ERROR")
		return -1
	end
	return poly_floor(poly_ceil(reverse(a)*c,n+m+1),m)
end
function t_kar_mul(n::Int64,a::Nemo.nmod_poly,c::Nemo.nmod_poly,err="")
	m = degree(a)
	err *= "n is $n m is $m deg(c) is $(degree(c))\n"
	err *= "c is $c\n"
	err *= "a is $a\n"
	if(degree(c)> m + n)
		println("ERROR in transpose mult")
		println(err)
		return zero(Z)
	end
	if n==0
		return reduce(+,zero(Z), [coeff(a,k)*coeff(c,k) for k in 0:m]) 
	end
	if m<=0
		return lead(a) * c
	end
	mu = div(m,2)+1
	vu = div(n,2)+1
	h = max(mu,vu)
	rho = h-1
	#rho = h
	if mu > n
		u = t_kar_mul(n,poly_ceil(a,mu),poly_ceil(c,n+rho),err)
		v = t_kar_mul(n,poly_floor(a,mu),poly_floor(c,mu),err)
		return u+v * one(Z)
	end
	if vu > m
		u = t_kar_mul(vu-1,a,poly_ceil(c,m+vu),err)
		v = t_kar_mul(n-vu,a,poly_floor(c,vu),err)
		return Z(u + v)
	end
	lambda = degree(poly_ceil(a,h) + poly_floor(a,h))
	r = t_kar_mul(h-1,poly_ceil(a,h),poly_ceil(c,rho+h)-poly_floor(poly_ceil(c,rho+2*h),h),err)
	s = t_kar_mul(n-h,poly_floor(a,h),poly_floor(c,2*h)-poly_floor(poly_ceil(c,m+n-h+1),h),err)
	t = t_kar_mul(h-1,poly_ceil(a,h) + poly_floor(a,h), poly_floor(poly_ceil(c,lambda + 2*h),h),err)
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
	#qi=coeff(one(Z),0)
	negi=-1
	for i=2:n
		qs[i]=qs[i-1]*q[i-1]
		#qi = q*qi
		vs[i]=(negi)*fs[i]*qs[i]
		us[i]=us[i-1]*q[i]*inv(1-q[i])
		ws[i]=vs[i]*inv(us[i])
		zs[i]=negi*us[i]*inv(qs[i])
		negi*=-1
	end
	#println(reduce(+,0,[x^i*us[i+1]for i in 0:n-1]))
	#println(reduce(+,0,[x^i*ws[i+1]for i in 0:n-1]))
	#println("ws is $ws")
	#println("us is $us")
	a = reduce(+,0,[x^i*us[i+1] for i in 0:n-1])
	c = reduce(+,0,[x^i*ws[i+1] for i in 0:n-1])
	F = t_mul(n,a,c)
	#G = t_kar_mul(n,a,c)
	#println("G is $G")
	#println("n is $n")
	#println("q is $qs")
	#println("G is $G")
	println("F is $F")
	r = [zs[i+1]*coeff(F,i) for i in 0:n-1]
	#println(r)
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
	#qi = 1
	for i=2:n
		qs[i] = qs[i-1]*q[i-1]
		#qi*=q
		us[i] = us[i-1]*(q[i]-1)
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
	#println(m)
	N=2*n-1
	p = convert(Int64,modulus(Z))
	q=2
	qs=Array{Nemo.nmod}(N)
	qs[1]=lead(one(Z))
	for qi = 2:p
		k = 2
		qk = lead(Z(qi))
		#qs[k]=qk
		while (qk != 1 && k <= N)
			qs[k] = qk
			qk=qk*qi
			k+=1
		end
		if k > N
			q = qi
			break
		end
	end
	#println("q is $qs")
	G = typeof(g)[]
	for i=1:m
		gi = g % x^(m) 
		g = divexact((g- gi),x^(m))
		push!(G, gi)
	end
	#println("G is $G")
	Λ = Array{Any}(m,N)
	H = Array{Any}(m)
	H[1] = one(Z)
	for i=2:m
		H[i]=(H[i-1]*h) % f
	end
	#println("H is $H")
	#A = fs_to_mat(H,n,m)
	A = Array{Any}(m,N)
	#qs are powers of q
	for i=1:m
		#qj=1
		for j=1:N
			A[i,j]=H[i](qs[j])
			#qj = q*qj
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
	#println("C is $BA")
	r = reduce(+, 0, [BA[i,:] .* Λ[i,:] for i in 1:m])
	#println("r is $r")
	coeffs =  newton_to_monomial_geom(newton_interp_geom(r,qs),qs)
	#println("coeffs $coeffs")
	rho = zero(Z)
	for i=0:length(coeffs)-1
		rho += x^i * coeffs[i+1]
	end
	#println("rho is $rho")
	return rho%f
end

#f = x^5+x^4+2x^2+3
#g = x^4+x^3+x+1
#h = x^4+32x^3+x^2+2
f = x^9+2x^3+2x+5
g = x^8+21x^3+66x+1
h = 3x^3+3x+9
#f = 14*x^13+13*x^12+12*x^11+11*x^10+10*x^9+9*x^8+8*x^7+7*x^6+6*x^5+5*x^4+4*x^3+3*x^2+2*x+1
#g = 93*x^12+94*x^11+27*x^10+25*x^9+43*x^8+84*x^7+26*x^6+62*x^5+52*x^4+51*x^3+22*x^2+102*x+13
#h = 16*x^11+15*x^10+x^9+72*x^8+23*x^7+77*x^6+73*x^5+13*x^4+27*x^3+15*x^2+26*x+15
T = 200
fails = 0
for i=0:T
	#println("f is $f")
	#println("g is $g")
	#println("h is $h")
	println("n is $(degree(f))")
	@time p = g(h)%f
	@time r = fast_comp_mod(f,g,h)
	if(r!=p)
		println("failed")
		println("fast $r")
		println()
		println("classical $p")
		fails++
		break
	end
	f+=x^(degree(f)+1)
	if(true)
		g+=x^(degree(g)+1)
	end
	h+=x^(degree(h)+1)
end
println("finished")
println("$(T-fails)/$T passed")
a = 2+3x+5x^2
c = 1+2x+3x^2+4x^3+5x^4
println(t_kar_mul(2,c,a))
println(t_mul(2,c,a))
