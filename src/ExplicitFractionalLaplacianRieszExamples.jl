module ExplicitFractionalLaplacianRieszExamples

using ClassicalOrthogonalPolynomials: hermiteh, jacobip, laguerrel
using HypergeometricFunctions: pFq, _₁F₁, _₂F₁
using MeijerG: meijerg
using QuadGK: quadgk
using SpecialFunctions: besselj, bessely, gamma

export PaperGFormula, evaluate, residue_reference, formula_1d, formula_highd,
       special_1d, special_highd, row5_input,
       riesz_row5_quadrature, fractional_laplacian_row1_exterior,
       logarithmic_potential_row1, logarithmic_potential_hermite,
       source_formula_1d, source_classical_1d, hypergeometric_1d,
       hypergeometric_highd, fractional_laplacian_source_exterior,
       higher_input_formula, higher_input_classical, higher_input_z_reduction,
       riesz_source_quadrature, fractional_laplacian_source_quadrature,
       higher_source_radial, riesz_highd_quadrature,
       fractional_laplacian_highd_separated_quadrature,
       fractional_laplacian_highd_quadrature

"""Meijer-G formula and prefactor."""
struct PaperGFormula{T}
    row::String
    prefactor::T
    a::Vector{T}
    b::Vector{T}
    m::Int
    n::Int
    z::T
end

evaluate(f::PaperGFormula) = f.prefactor * meijerg(f.a, f.b, f.m, f.n, f.z)

_reciprocal_gamma(x) = x <= 0 && isinteger(x) ? zero(x) : inv(gamma(x))

function _cancel_orders(a, b, m, n)
    al, ar = collect(a[1:n]), collect(a[(n+1):end])
    bl, br = collect(b[1:m]), collect(b[(m+1):end])
    changed = true
    while changed
        changed = false
        for i in eachindex(al)
            j = findfirst(isequal(al[i]), br)
            if j !== nothing
                deleteat!(al, i); deleteat!(br, j); n -= 1
                changed = true
                break
            end
        end
        changed && continue
        for i in eachindex(ar)
            j = findfirst(isequal(ar[i]), bl)
            if j !== nothing
                deleteat!(ar, i); deleteat!(bl, j); m -= 1
                changed = true
                break
            end
        end
    end
    vcat(al, ar), vcat(bl, br), m, n
end

function _lower_residues(a, b, m, n, z)
    p, q = length(a), length(b)
    argument = isodd(p - m - n) ? -z : z
    total = zero(complex(z))
    for k in 1:m
        bk = b[k]
        coefficient = z^bk
        for j in 1:m
            j == k || (coefficient *= gamma(b[j] - bk))
        end
        for j in 1:n
            coefficient *= gamma(1 + bk - a[j])
        end
        for j in (m + 1):q
            coefficient *= _reciprocal_gamma(1 + bk - b[j])
        end
        for j in (n + 1):p
            coefficient *= _reciprocal_gamma(a[j] - bk)
        end
        upper = Tuple(1 + bk - aj for aj in a)
        lower = Tuple(1 + bk - b[j] for j in 1:q if j != k)
        total += coefficient * pFq(upper, lower, argument)
    end
    total
end

"""Independent Slater-residue evaluation."""
function residue_reference(f::PaperGFormula)
    a, b, m, n = _cancel_orders(f.a, f.b, f.m, f.n)
    p, q = length(a), length(b)
    value = if p == q && f.z > one(f.z)
        aa = vcat(one(eltype(b)) .- b[1:m], one(eltype(b)) .- b[(m + 1):end])
        bb = vcat(one(eltype(a)) .- a[1:n], one(eltype(a)) .- a[(n + 1):end])
        _lower_residues(aa, bb, n, m, inv(f.z))
    else
        _lower_residues(a, b, m, n, f.z)
    end
    f.prefactor * value
end

_parity(n) = n - 2fld(n, 2)
_pochhammer(a, n) = gamma(a + n) / gamma(a)
_jacobi_at_one(alpha, n) = _pochhammer(alpha + 1, n) / factorial(n)

function row5_input(x::Real; a::Real, b::Real, n::Integer)
    abs(x) < 1 || return zero(float(x))
    (1-x^2)^a * (x^2)^b * jacobip(n, a, b, 2x^2-1)
end

function _row1_input(x::Real; a::Real, n::Integer)
    abs(x) < 1 || return zero(float(x))
    (1-x^2)^a * jacobip(n, a, a, x)
end

"""Table 5 Rows 1*--1** exterior check."""
function fractional_laplacian_row1_exterior(x::Real; s::Real, a::Real,
                                            n::Integer, rtol::Real=2e-10)
    s > 0 || throw(ArgumentError("s must be positive"))
    abs(x) > 1 || throw(ArgumentError("this helper requires |x| > 1"))
    integrand(y) = -_row1_input(y; a, n) / abs(x-y)^(1+2s)
    integral = quadgk(integrand, -1.0, 0.0, 1.0; rtol)[1]
    4^s * gamma(1/2+s) / (sqrt(pi)*abs(gamma(-s))) * integral
end

"""Table 5 Row 1*** logarithmic check."""
function logarithmic_potential_row1(x::Real; a::Real, n::Integer,
                                    rtol::Real=2e-10)
    n >= 1 || throw(ArgumentError("the logarithmic normalization requires n >= 1"))
    integrand(y) = _row1_input(y; a, n) * log(abs(x-y))
    points = -1 < x < 1 ? sort([-1.0, 0.0, float(x), 1.0]) : [-1.0, 0.0, 1.0]
    integral = sum(quadgk(integrand, points[i], points[i+1]; rtol)[1]
                   for i in 1:(length(points)-1) if points[i] != points[i+1])
    -integral/pi
end

"""Table 5 Row 6** logarithmic check."""
function logarithmic_potential_hermite(x::Real; n::Integer, rtol::Real=2e-10)
    n >= 1 || throw(ArgumentError("the logarithmic normalization requires n >= 1"))
    integrand(y) = exp(-y^2) * hermiteh(n, y) * log(abs(x-y))
    integral = quadgk(integrand, -Inf, x; rtol)[1] +
               quadgk(integrand, x, Inf; rtol)[1]
    -integral/pi
end

"""Table 5 Row 5 Riesz quadrature."""
function riesz_row5_quadrature(x::Real; t::Real, a::Real, b::Real, n::Integer,
                              rtol::Real=2e-10)
    0 < t < 1/2 || throw(ArgumentError("t must lie in (0, 1/2)"))
    -1 < x < 1 || throw(ArgumentError("this quadrature helper requires -1 < x < 1"))
    power = inv(2t)
    left(u) = row5_input(x-u^power; a, b, n) / (2t)
    right(u) = row5_input(x+u^power; a, b, n) / (2t)
    left_end, right_end = (x+1)^(2t), (1-x)^(2t)
    left_points = x > 0 ? [0.0, x^(2t), left_end] : [0.0, left_end]
    right_points = x < 0 ? [0.0, (-x)^(2t), right_end] : [0.0, right_end]
    integrate(g, points) = sum(quadgk(g, points[i], points[i+1]; rtol)[1]
                               for i in 1:(length(points)-1))
    integral = integrate(left, left_points) + integrate(right, right_points)
    gamma(1/2-t) / (4^t*sqrt(pi)*gamma(t)) * integral
end

"""Table 1 Meijer-G input."""
function source_formula_1d(row::Integer, x::Real; n::Integer=0,
                           a::Real=0.7, b::Real=0.4, lambda::Real=0.8,
                           alpha::Real=0.3, nu::Real=0.4, mu::Real=0.3)
    T = promote_type(Float64, typeof(float(x)), typeof(float(a)), typeof(float(b)),
                     typeof(float(lambda)), typeof(float(alpha)), typeof(float(nu)),
                     typeof(float(mu)))
    x, a, b, lambda, alpha, nu, mu = T.((x, a, b, lambda, alpha, nu, mu))
    k, ell, z = fld(n, 2), _parity(n), x^2
    if row in 1:4
        parameter = row == 1 ? a+1 : row == 2 ? lambda+T(1//2) :
                    row == 3 ? T(1//2) : T(3//2)
        aa = [parameter+k, -n+k+T(1//2)]
        bb = [zero(T), -n+2k+T(1//2)]
        c = if row == 1
            gamma(a+n+1)/factorial(n)
        elseif row == 2
            gamma(lambda+n+T(1//2))*_pochhammer(2lambda,n) /
            (factorial(n)*_pochhammer(lambda+T(1//2),n))
        elseif row == 3
            gamma(n+T(1//2))/(factorial(n)*_jacobi_at_one(T(-1//2),n))
        else
            (n+1)*gamma(n+T(3//2))/(factorial(n)*_jacobi_at_one(T(1//2),n))
        end
        return PaperGFormula(string(row), c*x^ell, aa, bb, 2, 0, z)
    elseif row == 5
        return PaperGFormula("5", gamma(a+n+1)/factorial(n),
                             [a+b+n+1, -T(n)], [zero(T), b], 2, 0, z)
    elseif row == 6
        return PaperGFormula("6", T(2)^n*x^ell,
                             [-n+k+T(1//2)], [zero(T), -n+2k+T(1//2)], 2, 0, z)
    elseif row == 7
        return PaperGFormula("7", inv(factorial(n)), [-n-alpha],
                             [zero(T), -alpha], 1, 1, z)
    elseif row == 8
        return PaperGFormula("8", one(T), T[], [nu/2, -nu/2], 1, 0, z)
    elseif row == 9 || row == 10
        phase = row == 9 ? a/T(pi)+(nu+1)/2 : a/T(pi)+nu/2
        return PaperGFormula(string(row), inv(sqrt(T(2))),
                             [T(1//4), T(3//4), phase],
                             [nu/2, (nu+1)/2, -nu/2, (1-nu)/2, phase], 2, 2, z)
    elseif row == 11
        return PaperGFormula("11", inv(sqrt(T(pi))), [zero(T), T(1//2)],
                             [(mu+nu)/2, -(mu+nu)/2, (mu-nu)/2, (nu-mu)/2],
                             1, 2, z)
    elseif row == 12
        return PaperGFormula("12", one(T), [-(nu+1)/2],
                             [nu/2, -nu/2, -(nu+1)/2], 2, 0, z)
    elseif row == 13
        return PaperGFormula("13", inv(sqrt(T(2))),
                             [T(1//4), T(3//4), -(nu+1)/2],
                             [-nu/2, nu/2, -(nu+1)/2, (1-nu)/2, (nu+1)/2],
                             2, 2, z)
    elseif row == 14
        return PaperGFormula("14", inv(sqrt(T(2))),
                             [T(1//4), T(3//4), -nu/2],
                             [(nu+1)/2, (1-nu)/2, -nu/2, -nu/2, nu/2],
                             2, 2, z)
    end
    throw(ArgumentError("row must be between 1 and 14"))
end

"""Table 1 classical input."""
function source_classical_1d(row::Integer, x::Real; n::Integer=0,
                             a::Real=0.7, b::Real=0.4, lambda::Real=0.8,
                             alpha::Real=0.3, nu::Real=0.4, mu::Real=0.3)
    ax = abs(x)
    if row == 1
        return ax < 1 ? (1-x^2)^a*jacobip(n,a,a,x) : zero(float(x))
    elseif row == 2
        gegenbauer = _pochhammer(2lambda,n)/_pochhammer(lambda+1/2,n) *
                     jacobip(n,lambda-1/2,lambda-1/2,x)
        return ax < 1 ? (1-x^2)^(lambda-1/2)*gegenbauer : zero(float(x))
    elseif row == 3
        return ax < 1 ? cos(n*acos(x))/sqrt(1-x^2) : zero(float(x))
    elseif row == 4
        return ax < 1 ? sqrt(1-x^2)*sin((n+1)*acos(x))/sqrt(1-x^2) : zero(float(x))
    elseif row == 5
        return row5_input(x; a, b, n)
    elseif row == 6
        return exp(-x^2)*hermiteh(n,x)
    elseif row == 7
        return exp(-x^2)*laguerrel(n,alpha,x^2)
    elseif row == 8
        return besselj(nu,2ax)
    elseif row == 9
        return cos(a+ax)*besselj(nu,ax)
    elseif row == 10
        return sin(a+ax)*besselj(nu,ax)
    elseif row == 11
        return besselj(mu,ax)*besselj(nu,ax)
    elseif row == 12
        return bessely(nu,2ax)
    elseif row == 13
        return cos(ax)*bessely(nu,ax)
    elseif row == 14
        return sin(ax)*bessely(nu,ax)
    end
    throw(ArgumentError("row must be between 1 and 14"))
end

"""Table 1 Rows 1--5 exterior quadrature."""
function fractional_laplacian_source_exterior(row::Integer, x::Real; s::Real,
                                              n::Integer=0, a::Real=0.7,
                                              b::Real=0.4, lambda::Real=0.8,
                                              rtol::Real=2e-10)
    row in 1:5 || throw(ArgumentError("row must be between 1 and 5"))
    s > 0 || throw(ArgumentError("s must be positive"))
    abs(x) > 1 || throw(ArgumentError("this helper requires |x| > 1"))
    integrand(y) = -source_classical_1d(row,y;n,a,b,lambda) / abs(x-y)^(1+2s)
    integral = quadgk(integrand,-1.0,0.0,1.0;rtol)[1]
    4^s*gamma(1/2+s)/(sqrt(pi)*abs(gamma(-s)))*integral
end

"""Table 1 fractional-Laplacian quadrature."""
function fractional_laplacian_source_quadrature(row::Integer, x::Real; s::Real,
                                                 n::Integer=0, a::Real=0.7,
                                                 b::Real=0.4, lambda::Real=0.8,
                                                 alpha::Real=0.3, nu::Real=0.4,
                                                 mu::Real=0.3,
                                                 rtol::Real=3e-8,
                                                 tail_cutoff::Real=256.0)
    row in 1:14 || throw(ArgumentError("row must be between 1 and 14"))
    0 < s < 1 || throw(ArgumentError("this quadrature helper requires 0 < s < 1"))
    x != 0 || throw(ArgumentError("use a nonzero evaluation point"))
    f(y) = source_classical_1d(row, y; n, a, b, lambda, alpha, nu, mu)
    fx = f(x)
    integrand(h) = (2fx - f(x-h) - f(x+h)) / h^(1+2s)

    # Split wherever either translated argument meets the origin or, for the
    # compactly supported rows, an endpoint of the support.
    breaks = Float64[0.0, abs(float(x))]
    if row <= 5
        append!(breaks, abs.((float(x)-1, float(x)+1)))
    end
    filter!(>(0.0), breaks)
    sort!(breaks)
    unique!(breaks)
    pushfirst!(breaks, 0.0)
    tail_start = max(breaks[end], row <= 5 ? breaks[end] : 8.0)
    tail_cutoff > tail_start || throw(ArgumentError("tail_cutoff is too small"))
    finite_breaks = sort!(unique!(vcat(breaks,tail_start)))
    total = sum(quadgk(integrand, finite_breaks[i], finite_breaks[i+1]; rtol)[1]
                for i in 1:(length(finite_breaks)-1)
                if finite_breaks[i] != finite_breaks[i+1])
    total += fx*tail_start^(-2s)/s
    if row > 5
        source_tail(h) = -(f(x-h)+f(x+h))/h^(1+2s)
        total += quadgk(source_tail,tail_start,tail_cutoff;rtol)[1]
    end
    4^s * gamma(1/2+s) / (sqrt(pi)*abs(gamma(-s))) * total
end

"""Table 1 Rows 1--7 Riesz quadrature."""
function riesz_source_quadrature(row::Integer, x::Real; t::Real,
                                 n::Integer=0, a::Real=0.7, b::Real=0.4,
                                 lambda::Real=0.8, alpha::Real=0.3,
                                 rtol::Real=3e-9)
    row in 1:7 || throw(ArgumentError("row must be between 1 and 7"))
    0 < t < 1/2 || throw(ArgumentError("t must lie in (0,1/2)"))
    power = inv(2t)
    left(u) = source_classical_1d(row,x-u^power;n,a,b,lambda,alpha)/(2t)
    right(u) = source_classical_1d(row,x+u^power;n,a,b,lambda,alpha)/(2t)
    if row <= 5 && -1 < x < 1
        left_points = x > 0 ? [0.0,x^(2t),(x+1)^(2t)] : [0.0,(x+1)^(2t)]
        right_points = x < 0 ? [0.0,(-x)^(2t),(1-x)^(2t)] : [0.0,(1-x)^(2t)]
        integrate(g,points) = sum(quadgk(g,points[i],points[i+1];rtol)[1]
                                  for i in 1:(length(points)-1))
        integral = integrate(left,left_points)+integrate(right,right_points)
    elseif row <= 5 && x > 1
        points = [(x-1)^(2t),x^(2t),(x+1)^(2t)]
        integral = sum(quadgk(left,points[i],points[i+1];rtol)[1]
                       for i in 1:(length(points)-1))
    elseif row <= 5 && x < -1
        points = [(-x-1)^(2t),(-x)^(2t),(-x+1)^(2t)]
        integral = sum(quadgk(right,points[i],points[i+1];rtol)[1]
                       for i in 1:(length(points)-1))
    else
        integral = quadgk(left,0.0,Inf;rtol)[1]+quadgk(right,0.0,Inf;rtol)[1]
    end
    gamma(1/2-t)/(4^t*sqrt(pi)*gamma(t))*integral
end

"""Table 6 radial input."""
function higher_source_radial(row::Symbol, r::Real; n::Integer=0,
                              a::Real=0.7, b::Real=0.3,
                              alpha::Real=0.2)
    if row == :A
        return r < 1 ? (1-r^2)^a*jacobip(n,a,b,2r^2-1) : 0.0
    elseif row == :B
        return exp(-r^2)*laguerrel(n,alpha,r^2)
    elseif row == :C
        return r > 1 ? (r^2-1)^a*jacobip(n,a,b,2r^2-1) : 0.0
    elseif row == :D
        return r > 1 ? (r^2-1)^a*jacobip(n,a,b,2/r^2-1) : 0.0
    end
    throw(ArgumentError("row must be :A, :B, :C, or :D"))
end

_sphere_area(d::Integer) = 2*pi^(d/2)/gamma(d/2)

function _hyp2f1_riesz_kernel(a::Real,b::Real,c::Real,z::Real,t::Real)
    delta = 2t-1 # c-a-b for the parameters used below
    (z < 0.8 || abs(delta-round(delta)) < 1e-12) && return _₂F₁(a,b,c,z)
    w = 1-z
    first = gamma(c)*gamma(delta)/(gamma(c-a)*gamma(c-b)) *
            _₂F₁(a,b,1-delta,w)
    second = w^delta * gamma(c)*gamma(-delta)/(gamma(a)*gamma(b)) *
             _₂F₁(c-a,c-b,1+delta,w)
    first+second
end

"""Table 6 angular kernel with `V_ell(e_1)=1`."""
function _riesz_angular_kernel(r::Real, rho::Real; t::Real, d::Integer,
                               ell::Integer)
    r > 0 || throw(ArgumentError("r must be positive"))
    rho >= 0 || throw(ArgumentError("rho must be nonnegative"))
    rho == 0 && return ell == 0 ? _sphere_area(d)*r^(-d+2t) : 0.0
    scale = max(r,rho)
    q = min(r,rho)/scale
    coefficient = _sphere_area(d) * _pochhammer(d/2-t,ell) /
                  _pochhammer(d/2,ell)
    aa, bb, cc = d/2-t+ell, 1-t, d/2+ell
    coefficient * q^ell * scale^(-d+2t) *
    _hyp2f1_riesz_kernel(aa,bb,cc,q^2,t)
end

"""Table 6 Riesz quadrature with `V_ell(r e_1)=r^ell`."""
function riesz_highd_quadrature(row::Symbol, r::Real; t::Real, d::Integer,
                                ell::Integer=0, n::Integer=0,
                                a::Real=0.7, b::Real=0.3,
                                alpha::Real=0.2, rtol::Real=2e-7,
                                tail_cutoff::Real=4096.0)
    row in (:A,:B,:C,:D) || throw(ArgumentError("invalid Table 6 row"))
    0 < t < d/2 || throw(ArgumentError("t must lie in (0,d/2)"))
    r > 0 || throw(ArgumentError("r must be positive"))
    radial(rho) = rho^(d-1+ell) *
                  higher_source_radial(row,rho;n,a,b,alpha) *
                  _riesz_angular_kernel(r,rho;t,d,ell)
    lower, upper = row == :A ? (0.0,1.0) :
                   row == :B ? (0.0,min(12.0,float(tail_cutoff))) :
                   (1.0,float(tail_cutoff))
    integral = if row in (:C,:D) && r < 1
        # Remove the (rho-1)^a endpoint singularity of the exterior inputs.
        power = inv(a+1)
        function boundary(u)
            delta = u^power
            rho = 1+delta
            argument = row == :C ? 2rho^2-1 : 2/rho^2-1
            source = (delta*(rho+1))^a*jacobip(n,a,b,argument)
            rho^(d-1+ell)*source*_riesz_angular_kernel(r,rho;t,d,ell)*
            power*u^(power-1)
        end
        quadgk(boundary,0.0,1.0;rtol)[1] +
        quadgk(radial,2.0,upper;rtol)[1]
    elseif lower < r < upper
        # At rho=r the angular integral behaves like |rho-r|^(2t-1).
        # The following substitutions remove that integrable singularity.
        power = inv(2t)
        left(u) = radial(r-u^power)*power*u^(power-1)
        right(u) = radial(r+u^power)*power*u^(power-1)
        quadgk(left,0.0,(r-lower)^(2t);rtol)[1] +
        quadgk(right,0.0,(upper-r)^(2t);rtol)[1]
    else
        quadgk(radial,lower,upper;rtol)[1]
    end
    gamma(d/2-t)/(4^t*pi^(d/2)*gamma(t)) * integral
end

"""Table 6 separated-support quadrature."""
function fractional_laplacian_highd_separated_quadrature(
        row::Symbol, r::Real; s::Real, d::Integer, ell::Integer=0,
        n::Integer=0, a::Real=0.7, b::Real=0.3,
        alpha::Real=0.2, rtol::Real=2e-7,
        tail_cutoff::Real=4096.0)
    row in (:A,:C,:D) || throw(ArgumentError("row must be :A, :C, or :D"))
    0 < s < 1 || throw(ArgumentError("this helper requires 0 < s < 1"))
    separated = row == :A ? r > 1 : 0 < r < 1
    separated || throw(ArgumentError("evaluation point must be separated from the support"))
    radial(rho) = rho^(d-1+ell) *
                  higher_source_radial(row,rho;n,a,b,alpha) *
                  _riesz_angular_kernel(r,rho;t=-s,d,ell)
    points = row == :A ? [0.0,1.0] : [1.0,float(tail_cutoff)]
    integral = if row in (:C,:D)
        power = inv(a+1)
        function boundary(u)
            delta = u^power
            rho = 1+delta
            argument = row == :C ? 2rho^2-1 : 2/rho^2-1
            source = (delta*(rho+1))^a*jacobip(n,a,b,argument)
            rho^(d-1+ell)*source*_riesz_angular_kernel(r,rho;t=-s,d,ell)*
            power*u^(power-1)
        end
        quadgk(boundary,0.0,1.0;rtol)[1] +
        quadgk(radial,2.0,tail_cutoff;rtol)[1]
    else
        sum(quadgk(radial,points[i],points[i+1];rtol)[1]
            for i in 1:(length(points)-1))
    end
    -4^s*gamma(d/2+s)/(pi^(d/2)*abs(gamma(-s))) * integral
end

function _zonal_solid_harmonic(ell::Integer, d::Integer, rho::Real,
                               cosine::Real)
    ell == 0 && return one(float(rho))
    rho == 0 && return zero(float(rho))
    c = clamp(cosine, -one(float(cosine)), one(float(cosine)))
    angular = if d == 2
        cos(ell*acos(c))
    else
        alpha = (d-3)/2
        jacobip(ell,alpha,alpha,c)/_jacobi_at_one(alpha,ell)
    end
    rho^ell*angular
end

"""Table 6 fractional-Laplacian quadrature with `V_ell(r e_1)=r^ell`."""
function fractional_laplacian_highd_quadrature(
        row::Symbol, r::Real; s::Real, d::Integer, ell::Integer=0,
        n::Integer=0, a::Real=0.7, b::Real=0.3,
        alpha::Real=0.2, rtol::Real=2e-5,
        tail_cutoff::Real=32.0)
    row in (:A,:B,:C,:D) || throw(ArgumentError("invalid Table 6 row"))
    0 < s < 1 || throw(ArgumentError("this helper requires 0 < s < 1"))
    d >= 2 || throw(ArgumentError("this helper requires d >= 2"))
    ell >= 0 || throw(ArgumentError("ell must be nonnegative"))
    r > 0 || throw(ArgumentError("r must be positive"))

    source(rho, cosine) = _zonal_solid_harmonic(ell,d,rho,cosine) *
                          higher_source_radial(row,rho;n,a,b,alpha)
    fx = r^ell*higher_source_radial(row,r;n,a,b,alpha)
    angular_rtol = max(rtol/5,1e-8)

    function angular_points(h)
        points = Float64[0.0,pi]
        if row != :B && h > 0
            u = (1-r^2-h^2)/(2r*h)
            if -1 < u < 1
                push!(points,acos(u),acos(-u))
            end
        end
        sort!(points)
        unique!(points)
    end

    function angular_second_difference(h)
        h == 0 && return zero(float(h))
        function integrand(theta)
            u = cos(theta)
            rho_plus = sqrt(max(0.0,r^2+h^2+2r*h*u))
            rho_minus = sqrt(max(0.0,r^2+h^2-2r*h*u))
            cosine_plus = rho_plus == 0 ? 1.0 : (r+h*u)/rho_plus
            cosine_minus = rho_minus == 0 ? 1.0 : (r-h*u)/rho_minus
            (2fx-source(rho_plus,cosine_plus)-source(rho_minus,cosine_minus)) *
            sin(theta)^(d-2)
        end
        _sphere_area(d-1)*quadgk(integrand,angular_points(h)...;
                                 rtol=angular_rtol)[1] /
        h^(1+2s)
    end

    function angular_source_tail(h)
        function integrand(theta)
            u = cos(theta)
            rho_plus = sqrt(max(0.0,r^2+h^2+2r*h*u))
            rho_minus = sqrt(max(0.0,r^2+h^2-2r*h*u))
            cosine_plus = rho_plus == 0 ? 1.0 : (r+h*u)/rho_plus
            cosine_minus = rho_minus == 0 ? 1.0 : (r-h*u)/rho_minus
            -(source(rho_plus,cosine_plus)+source(rho_minus,cosine_minus)) *
            sin(theta)^(d-2)
        end
        _sphere_area(d-1)*quadgk(integrand,angular_points(h)...;
                                 rtol=angular_rtol)[1] /
        h^(1+2s)
    end

    support_end = row == :A ? r+1 : max(r+1,8.0)
    tail_cutoff > support_end || throw(ArgumentError("tail_cutoff is too small"))
    breaks = Float64[0.0,r,abs(r-1),r+1,support_end]
    filter!(x -> 0 <= x <= support_end,breaks)
    sort!(breaks)
    unique!(breaks)
    total = sum(quadgk(angular_second_difference,breaks[i],breaks[i+1];rtol)[1]
                for i in 1:(length(breaks)-1)
                if breaks[i] != breaks[i+1])

    # Integrate the constant part of the infinite tail exactly.  Row A is
    # compactly supported; the other source tails are truncated explicitly.
    total += _sphere_area(d)*fx*support_end^(-2s)/s
    if row != :A
        total += quadgk(angular_source_tail,support_end,tail_cutoff;rtol)[1]
    end
    4^s*gamma(d/2+s)/(2pi^(d/2)*abs(gamma(-s))) * total
end

"""Jacobi inputs preceding Table 6."""
function higher_input_formula(row::Symbol, r::Real; n::Integer=0,
                              a::Real=0.7, b::Real=0.3, D::Real=1.5)
    z, c = r^2, gamma(a+n+1)/factorial(n)
    if row == :A
        return PaperGFormula("input A",c,[a+n+1,-n-b],[-b,0.0],2,0,z)
    elseif row == :Z
        return PaperGFormula("input Z",(-1)^n*c,[1-D-n,a+n+1],[0.0,1-D],1,1,z)
    elseif row == :C
        return PaperGFormula("input C",c,[-b-n,a+n+1],[0.0,-b],0,2,z)
    elseif row == :D
        return PaperGFormula("input D",c,[a+1,a+b+1],[-float(n),a+b+n+1],0,2,z)
    end
    throw(ArgumentError("row must be :A, :Z, :C, or :D"))
end

function higher_input_classical(row::Symbol, r::Real; n::Integer=0,
                                a::Real=0.7, b::Real=0.3, D::Real=1.5)
    if row == :A || row == :Z
        beta = row == :Z ? D-1 : b
        return r < 1 ? (1-r^2)^a*jacobip(n,a,beta,2r^2-1) : 0.0
    elseif row == :C
        return r > 1 ? (r^2-1)^a*jacobip(n,a,b,2r^2-1) : 0.0
    elseif row == :D
        return r > 1 ? (r^2-1)^a*jacobip(n,a,b,2/r^2-1) : 0.0
    end
    throw(ArgumentError("row must be :A, :Z, :C, or :D"))
end

"""Hypergeometric reduction preceding Table 6."""
function higher_input_z_reduction(r::Real; n::Integer=0, a::Real=0.7,
                                  D::Real=1.5)
    (-1)^n*gamma(D+n)/(factorial(n)*gamma(D)) *
    _₂F₁(D+n,-a-n,D,r^2)
end

"""Table 5 special cases."""
function special_1d(row::AbstractString, x::Real; n::Integer=0,
                    a::Real=0.7, s::Real=0.3)
    k, ell, km = fld(n, 2), _parity(n), fld(n-1, 2)
    z, ax = x^2, abs(x)
    if row == "1*"
        if ax < 1
            return 4^s * gamma(s+k+1) * gamma(n+s-k+1/2) /
                   (factorial(k) * gamma(n-k+1/2)) * jacobip(n, s, s, x)
        end
        return -sin(pi*s) * x^ell * gamma(n+s+1) * gamma(n+2s+1) *
               ax^(-2km-2s-3) * _₂F₁(s+km+3/2, s+k+1,
                                         n+s+3/2, inv(z)) /
               (2^n * sqrt(pi) * factorial(n) * gamma(n+s+3/2))
    elseif row == "1**"
        common = gamma(a+n+1)/factorial(n) * x^ell
        if ax < 1
            value = 2*(-1)^k * factorial(fld(n+1, 2)) *
                    _₂F₁(-a-k+1/2, n-k+1, ell+1/2, z) /
                    (gamma(ell+1/2) * gamma(a+k+1/2))
            return common * value
        end
        value = -2.0^(-n) * gamma(n+2) * ax^(-2km-4) *
                _₂F₁(k+3/2, km+2, (2n+3)/2+a, inv(z)) /
                (sqrt(pi) * gamma((2n+3)/2+a))
        return common * value
    elseif row == "1***"
        n >= 1 || throw(ArgumentError("Row 1*** requires n >= 1"))
        if ax < 1
            return (-1)^k * gamma(a+n+1) * factorial(km) * x^ell *
                   _₂F₁(-a-k-1/2, n-k, ell+1/2, z) /
                   (2 * factorial(n) * gamma(ell+1/2) * gamma(a+k+3/2))
        end
        return 2.0^(-n) * gamma(n) * gamma(a+n+1) * x^ell *
               ax^(-2*(km+1)) * _₂F₁(km+1, k+1/2,
                                         a+n+3/2, inv(z)) /
               (sqrt(pi) * factorial(n) * gamma(a+n+3/2))
    elseif row == "6*"
        return 2.0^n * (-1)^k * 2.0^(-2k+n+1) * x^ell *
               gamma(fld(n+1, 2)+1) / sqrt(pi) *
               _₁F₁(fld(n+1, 2)+1, ell+1/2, -z)
    elseif row == "6**"
        n >= 1 || throw(ArgumentError("Row 6** requires n >= 1"))
        return (-1)^k * 2.0^(n-1) * factorial(km) * x^ell /
               gamma(ell+1/2) * _₁F₁(n-k, ell+1/2, -z)
    end
    throw(ArgumentError("row must be \"1*\", \"1**\", \"1***\", \"6*\", or \"6**\""))
end

"""Table 7 starred cases."""
function special_highd(row::AbstractString, r::Real; s::Real, d::Integer,
                       ell::Integer=0, n::Integer=0, a::Real=0.7,
                       solid_harmonic::Real=1.0)
    D, z, v = (d+2ell)/2, r^2, solid_harmonic
    if row == "A*"
        common = v * (-1)^n * 4^s * gamma(1+a+n) / factorial(n)
        if r < 1
            return common * gamma(D+n+s) *
                   _₂F₁(-a-n+s, D+n+s, D, z) /
                   (gamma(D) * gamma(a+n-s+1))
        end
        return common * gamma(D+n+s) * r^(-d-2ell-2n-2s) *
               _₂F₁(n+s+1, D+n+s, a+D+2n+1, inv(z)) /
               (gamma(-n-s) * gamma(a+D+2n+1))
    elseif row == "A**"
        common = v * 4^s * gamma(1+s+n) / factorial(n)
        if r < 1
            return common * gamma(D+n+s) /
                   gamma(d/2+n+ell) *
                   jacobip(n, s, D-1, 2z-1)
        end
        return common * (-1)^n * gamma(D+n+s) * r^(-d-2ell-2n-2s) *
               _₂F₁(n+s+1, D+n+s, s+D+2n+1, inv(z)) /
               (gamma(-n-s) * gamma(s+D+2n+1))
    elseif row == "A***"
        common = v * 4^s
        if r < 1
            return common * gamma(D+s) * _₂F₁(s, D+s, D, z) /
                   (gamma(D) * gamma(1-s))
        end
        return common * gamma(D+s) * r^(-d-2s-2ell) *
               _₂F₁(s+1, D+s, D+1, inv(z)) /
               (gamma(D+1) * gamma(-s))
    end
    throw(ArgumentError("row must be \"A*\", \"A**\", or \"A***\""))
end

"""Table 7 Rows A--B."""
function hypergeometric_highd(row::Symbol, r::Real; s::Real, d::Integer,
                              ell::Integer=0, n::Integer=0, a::Real=0.7,
                              b::Real=0.3, alpha::Real=0.2,
                              solid_harmonic::Real=1.0)
    D, z, v = (d+2ell)/2, r^2, solid_harmonic
    if row == :A
        common = v * 4^s * gamma(a+n+1)/factorial(n)
        if r < 1
            return common * gamma(-b-s)*gamma(D+s) *
                   pFq((-a-n+s,b+n+s+1,D+s),(b+s+1,D),z) /
                   (gamma(D)*gamma(a+n-s+1)*gamma(-b-n-s))
        end
        return common * gamma(-b+D)*gamma(D+s)*r^(-d-2*(s+ell)) *
               pFq((s+1,-b+D,D+s),(-b+D-n,a+D+n+1),inv(z)) /
               (gamma(-s)*gamma(a+D+n+1)*gamma(-b+D-n))
    elseif row == :B
        return v * 4^s * gamma(D+s)*gamma(n+s+alpha+1) /
               (factorial(n)*gamma(D)*gamma(s+alpha+1)) *
               pFq((D+s,n+s+alpha+1),(D,s+alpha+1),-z)
    end
    throw(ArgumentError("row must be :A or :B"))
end

"""Table 2 Meijer-G formulas."""
function formula_1d(row::Integer, x::Real; s::Real, n::Integer=0,
                    a::Real=0.7, b::Real=0.4, lambda::Real=0.8,
                    alpha::Real=0.3, nu::Real=0.4, mu::Real=0.3)
    T = promote_type(Float64, typeof(float(x)), typeof(float(s)), typeof(float(a)),
                     typeof(float(b)), typeof(float(lambda)), typeof(float(alpha)),
                     typeof(float(nu)), typeof(float(mu)))
    x, s, a, b, lambda, alpha, nu, mu = T.((x, s, a, b, lambda, alpha, nu, mu))
    k, ell, z = fld(n, 2), _parity(n), x^2
    if row == 1
        aa = [T(1//2)-s-ell, a+k+1-s, -n+k+T(1//2)-s]
        bb = [zero(T), -n+2k+T(1//2)-s, T(1//2)-ell]
        c = 4^s * gamma(a+n+1) / factorial(n) * x^ell
        return PaperGFormula("1", c, aa, bb, 2, 1, z)
    elseif row == 2
        aa = [T(1//2)-s-ell, lambda+T(1//2)+k-s, -n+k+T(1//2)-s]
        bb = [zero(T), -n+2k+T(1//2)-s, T(1//2)-ell]
        c = 4^s * gamma(lambda+T(1//2)+n) * _pochhammer(2lambda, n) /
            (factorial(n) * _pochhammer(lambda+T(1//2), n)) * x^ell
        return PaperGFormula("2", c, aa, bb, 2, 1, z)
    elseif row == 3
        aa = [T(1//2)-s-ell, k+T(1//2)-s, -n+k+T(1//2)-s]
        bb = [zero(T), -n+2k+T(1//2)-s, T(1//2)-ell]
        c = 4^s * gamma(n+T(1//2)) /
            (factorial(n) * _jacobi_at_one(T(-1//2), n)) * x^ell
        return PaperGFormula("3", c, aa, bb, 2, 1, z)
    elseif row == 4
        aa = [T(1//2)-s-ell, k+T(3//2)-s, -n+k+T(1//2)-s]
        bb = [zero(T), -n+2k+T(1//2)-s, T(1//2)-ell]
        c = 4^s * (n+1) * gamma(n+T(3//2)) /
            (factorial(n) * _jacobi_at_one(T(1//2), n)) * x^ell
        return PaperGFormula("4", c, aa, bb, 2, 1, z)
    elseif row == 5
        aa = [T(1//2)-s, a+b+n-s+1, -n-s]
        bb = [zero(T), b-s, T(1//2)]
        c = 4^s * gamma(a+n+1) / factorial(n)
        return PaperGFormula("5", c, aa, bb, 2, 1, z)
    elseif row == 6
        aa = [-n-s+2k+T(1//2), -n-s+k+T(1//2)]
        bb = [zero(T), -n-s+2k+T(1//2), -n+2k+T(1//2)]
        c = 2^n * 4^s * x^ell
        return PaperGFormula("6", c, aa, bb, 2, 1, z)
    elseif row == 7
        aa = [T(1//2)-s, -n-s-alpha]
        bb = [zero(T), -alpha-s, T(1//2)]
        return PaperGFormula("7", 4^s/factorial(n), aa, bb, 1, 2, z)
    elseif row == 8
        aa = [T(1//2)-s, -s]
        bb = [zero(T), nu/2-s, T(1//2), -nu/2-s]
        return PaperGFormula("8", 4^s, aa, bb, 2, 1, z)
    elseif row == 9 || row == 10
        phase = row == 9 ? a/T(pi)+(nu+1)/2-s : a/T(pi)+nu/2-s
        aa = [T(1//4)-s, T(1//2)-s, T(3//4)-s, -s, phase]
        bb = [zero(T), nu/2-s, (nu+1)/2-s, T(1//2),
              (1-nu)/2-s, -s-nu/2, phase]
        return PaperGFormula(string(row), 2^(2s-T(1//2)), aa, bb, 3, 3, z)
    elseif row == 11
        aa = [T(1//2)-s, -s, T(1//2)-s, -s]
        bb = [zero(T), (mu+nu)/2-s, -(mu+nu)/2-s,
              (mu-nu)/2-s, (nu-mu)/2-s, T(1//2)]
        return PaperGFormula("11", 4^s/sqrt(T(pi)), aa, bb, 2, 3, z)
    elseif row == 12
        aa = [T(1//2)-s, -(nu+1)/2-s, -s]
        bb = [zero(T), nu/2-s, -nu/2-s, -(nu+1)/2-s, T(1//2)]
        return PaperGFormula("12", 4^s, aa, bb, 3, 1, z)
    elseif row == 13
        aa = [T(1//2)-s, T(1//4)-s, T(3//4)-s, -(nu+1)/2-s, -s]
        bb = [zero(T), -nu/2-s, nu/2-s, -(nu+1)/2-s,
              (1-nu)/2-s, (nu+1)/2-s, T(1//2)]
        return PaperGFormula("13", 2^(2s-T(1//2)), aa, bb, 3, 3, z)
    elseif row == 14
        aa = [T(1//2)-s, T(1//4)-s, T(3//4)-s, -s-nu/2, -s]
        bb = [zero(T), (1-nu)/2-s, (1+nu)/2-s, -nu/2-s,
              -nu/2-s, nu/2-s, T(1//2)]
        return PaperGFormula("14", 2^(2s-T(1//2)), aa, bb, 3, 3, z)
    end
    throw(ArgumentError("row must be between 1 and 14"))
end

"""Tables 3--4 hypergeometric formulas."""
function hypergeometric_1d(row::Integer, x::Real; s::Real, n::Integer=0,
                           a::Real=0.7, b::Real=0.4, lambda::Real=0.8,
                           alpha::Real=0.3, nu::Real=0.4, mu::Real=0.3)
    k, ell, km = fld(n,2), _parity(n), fld(n-1,2)
    z, ax = x^2, abs(x)
    if row in 1:4
        sigma = row == 1 ? a : row == 2 ? lambda-1/2 : row == 3 ? -1/2 : 1/2
        common = if row == 1
            4^s*gamma(a+n+1)/factorial(n)*x^ell
        elseif row == 2
            4^s*gamma(lambda+1/2+n)*_pochhammer(2lambda,n) /
            (factorial(n)*_pochhammer(lambda+1/2,n))*x^ell
        elseif row == 3
            4^s*gamma(n+1/2)/(factorial(n)*_jacobi_at_one(-1/2,n))*x^ell
        else
            4^s*(n+1)*gamma(n+3/2)/(factorial(n)*_jacobi_at_one(1/2,n))*x^ell
        end
        if ax < 1
            branch = pi * _₂F₁(-sigma+s-k, n+s-k+1/2, ell+1/2, z) /
                     (sin(pi*(2k-n-s+1/2))*gamma(ell+1/2)*
                      gamma(-n-s+k+1/2)*gamma(sigma-s+k+1))
            return common*branch
        end
        outer_c = row == 1 ? (2n+3)/2+a : row == 2 ? n+1+lambda :
                  row == 3 ? n+1 : n+2
        branch = -2.0^(-n-2s)*sin(pi*s)*gamma(n+2s+1)*
                 ax^(-2km-2s-3) *
                 _₂F₁(s+k+1, km+3/2+s, outer_c, inv(z)) /
                 (sqrt(pi)*gamma(outer_c))
        return common*branch
    elseif row == 5
        common = gamma(a+n+1)/factorial(n)
        if ax < 1
            first = gamma(s+1/2)*gamma(b-s) *
                    pFq((s+1/2,-a-b-n+s,n+s+1),(1/2,-b+s+1),z) /
                    (sqrt(pi)*gamma(-n-s)*gamma(a+b+n-s+1))
            second = z^(b-s)*gamma(s-b)*gamma(b+1/2) *
                     pFq((b+1/2,-a-n,n+b+1),(b-s+1,b-s+1/2),z) /
                     (gamma(b-s+1/2)*gamma(a+n+1)*gamma(-b-n))
            return common*4^s*(first+second)
        end
        return common * -gamma(b+1/2)*sin(pi*s)*gamma(2s+1)*ax^(-2s-1) *
               pFq((b+1/2,s+1/2,s+1),(1/2-n,a+b+n+3/2),inv(z)) /
               (sqrt(pi)*gamma(1/2-n)*gamma(a+b+n+3/2))
    elseif row == 6
        return 2.0^n*pi*4^s*x^ell *
               _₁F₁(n+s-k+1/2,ell+1/2,-z) /
               (sin(pi/2*(4k-2n-2s+1))*gamma(ell+1/2)*
                gamma(-n-s+k+1/2))
    elseif row == 7
        return 4^s*gamma(s+1/2)*gamma(n+s+alpha+1) *
               pFq((s+1/2,n+s+alpha+1),(1/2,s+alpha+1),-z) /
               (sqrt(pi)*factorial(n)*gamma(s+alpha+1))
    elseif row == 8
        first = sin(pi*nu/2)*z^((nu-2s)/2) *
                pFq((nu/2+1/2,nu/2+1),
                    (-s+nu/2+1/2,-s+nu/2+1,nu+1),-z) /
                (sin(pi*(nu-2s)/2)*gamma(-2s+nu+1))
        second = sin(pi*s)*gamma(2s+1) *
                 pFq((s+1/2,s+1),(1/2,s-nu/2+1,s+nu/2+1),-z) /
                 (sin(pi*(nu-2s)/2)*gamma(s-nu/2+1)*gamma(s+nu/2+1))
        return first-second
    elseif row == 9 || row == 10
        fodd = pFq((nu/2+3/4,nu/2+1,nu/2+5/4,nu/2+3/2),
                   (3/2,-s+nu/2+1,-s+nu/2+3/2,nu+1,nu+3/2),-z)
        feven = pFq((nu/2+1/4,nu/2+1/2,nu/2+3/4,nu/2+1),
                    (1/2,-s+nu/2+1/2,-s+nu/2+1,nu+1/2,nu+1),-z)
        fs = pFq((s+1/4,s+1/2,s+3/4,s+1),
                 (1/2,s-nu/2+1/2,s-nu/2+1,s+nu/2+1/2,s+nu/2+1),-z)
        if row == 9
            bracket = (nu+1)*sin(a)*cos(pi*nu/2)*ax*fodd /
                      ((-nu+2s-1)*cos(pi/2*(nu-2s+2))) +
                      cos(a)*sin(pi*nu/2)*feven/sin(pi*(s-nu/2))
            tail = 4.0^(-s)*sin(pi*s)*gamma(4s+1)*
                   cos(a+pi*nu/2-pi*s)*fs /
                   (sin(pi*(s-nu/2))*cos(pi*(s-nu/2))*
                    gamma(2s-nu+1)*gamma(2s+nu+1))
        else
            bracket = (nu+1)*cos(a)*cos(pi*nu/2)*ax*fodd /
                      ((nu-2s+1)*cos(pi/2*(nu-2s+2))) +
                      sin(a)*sin(pi*nu/2)*feven/sin(pi*(s-nu/2))
            tail = 2.0^(1-2s)*sin(pi*s)*gamma(4s+1)*
                   sin(a+pi*nu/2-pi*s)*fs /
                   (sin(2pi*s-pi*nu)*gamma(2s-nu+1)*gamma(2s+nu+1))
        end
        return -2.0^(-nu)*ax^(nu-2s)/gamma(-2s+nu+1)*bracket + tail
    elseif row == 11
        fsum = pFq(((mu+nu)/2+1/2,(mu+nu)/2+1/2,(mu+nu)/2+1,(mu+nu)/2+1),
                   (mu+1,-s+(mu+nu)/2+1/2,-s+(mu+nu)/2+1,nu+1,mu+nu+1),-z)
        fs = pFq((s+1/2,s+1/2,s+1,s+1),
                 (1/2,s-(mu+nu)/2+1,s+(mu-nu)/2+1,
                  s+(nu-mu)/2+1,s+(mu+nu)/2+1),-z)
        first = -2.0^(-mu-nu)*sin(pi*(mu+nu)/2)*gamma(mu+nu+1)*
                ax^(mu+nu-2s)*fsum /
                (gamma(mu+1)*gamma(nu+1)*cos(pi*(mu+nu-2s+1)/2)*
                 gamma(-2s+mu+nu+1))
        second = -4.0^(-s)*sin(pi*s)*gamma(2s+1)^2*fs /
                 (sin(pi*(mu+nu-2s)/2)*gamma(s-(mu+nu)/2+1)*
                  gamma(s+(mu-nu)/2+1)*gamma(s+(nu-mu)/2+1)*
                  gamma(s+(mu+nu)/2+1))
        return first+second
    elseif row == 12
        fs = pFq((s+1/2,s+1),(1/2,s-nu/2+1,s+nu/2+1),-z)
        fminus = pFq((1/2-nu/2,1-nu/2),
                     (1-nu,-s-nu/2+1/2,-s-nu/2+1),-z)
        fplus = pFq((nu/2+1/2,nu/2+1),
                    (-s+nu/2+1/2,-s+nu/2+1,nu+1),-z)
        first = sqrt(pi)*2.0^(nu+2s+1)*sin(pi*s)*gamma(2s+1)*fs /
                (sin(pi*(nu-2s)/2)*sin(pi*nu/2+pi*s)*
                 gamma(-s-nu/2-1/2)*gamma(s-nu/2+1)*gamma(2s+nu+2))
        second = -ax^(-nu-2s)*fminus /
                 (2cos(pi*nu/2)*sin(pi*nu/2+pi*s)*gamma(-2s-nu+1))
        third = -sin(pi*nu/2)*cos(pi*nu)*ax^(nu-2s)*fplus /
                (sin(pi*nu)*sin(pi*(s-nu/2))*gamma(-2s+nu+1))
        return first+second+third
    elseif row == 13
        fminus = pFq((1/4-nu/2,1/2-nu/2,3/4-nu/2,1-nu/2),
                     (1/2,1/2-nu,1-nu,-s-nu/2+1/2,-s-nu/2+1),-z)
        fs = pFq((s+1/4,s+1/2,s+3/4,s+1),
                 (1/2,s-nu/2+1/2,s-nu/2+1,s+nu/2+1/2,s+nu/2+1),-z)
        fplus = pFq((nu/2+1/4,nu/2+1/2,nu/2+3/4,nu/2+1),
                    (1/2,-s+nu/2+1/2,-s+nu/2+1,nu+1/2,nu+1),-z)
        first = -2.0^(nu-1)*ax^(-nu-2s)*fminus /
                (cos(pi*nu/2)*sin(pi*nu/2+pi*s)*gamma(-2s-nu+1))
        second = -4.0^(-s)*sin(pi*s)*gamma(4s+1)*
                 cos(pi*nu/2+pi*s)*fs /
                 (sin(pi*(nu-2s)/2)*sin(pi*nu/2+pi*s)*
                  gamma(2s-nu+1)*gamma(2s+nu+1))
        third = -2.0^(-nu)*sin(pi*nu/2)*cos(pi*nu)*ax^(nu-2s)*fplus /
                (sin(pi*nu)*sin(pi*(s-nu/2))*gamma(-2s+nu+1))
        return first+second+third
    elseif row == 14
        fplus = pFq((nu/2+3/4,nu/2+1,nu/2+5/4,nu/2+3/2),
                    (3/2,-s+nu/2+1,-s+nu/2+3/2,nu+1,nu+3/2),-z)
        fminus = pFq((3/4-nu/2,1-nu/2,5/4-nu/2,3/2-nu/2),
                     (3/2,1-nu,3/2-nu,-s-nu/2+1,-s-nu/2+3/2),-z)
        fs = pFq((s+1/4,s+1/2,s+3/4,s+1),
                 (1/2,s-nu/2+1/2,s-nu/2+1,s+nu/2+1/2,s+nu/2+1),-z)
        first = 2.0^(-nu)*(nu+1)*cos(pi*nu/2)*cos(pi*nu)*
                ax^(nu-2s+1)*fplus /
                (sin(pi*nu)*cos(pi*(s-nu/2))*gamma(-2s+nu+2))
        second = -2.0^(nu-1)*(nu-1)*ax^(-nu-2s+1)*fminus /
                 (sin(pi*nu/2)*sin(pi*(nu+2s-1)/2)*gamma(-2s-nu+2))
        third = 4.0^(-s)*sin(pi*s)*gamma(4s+1)*sin(pi*nu/2+pi*s)*fs /
                (cos(pi*nu/2+pi*s)*cos(pi*(s-nu/2))*
                 gamma(2s-nu+1)*gamma(2s+nu+1))
        return first+second+third
    end
    throw(ArgumentError("row must be between 1 and 14"))
end

"""Table 6 Meijer-G formulas."""
function formula_highd(row::Symbol, r::Real; s::Real, d::Integer, ell::Integer=0,
                       n::Integer=0, a::Real=0.7, b::Real=0.3,
                       alpha::Real=0.2, solid_harmonic::Real=1.0)
    T = promote_type(Float64, typeof(float(r)), typeof(float(s)), typeof(float(a)),
                     typeof(float(b)), typeof(float(alpha)), typeof(float(solid_harmonic)))
    r, s, a, b, alpha, v = T.((r, s, a, b, alpha, solid_harmonic))
    D, z = T(d+2ell)/2, r^2
    c = v * 4^s / factorial(n)
    if row == :A
        aa = [1-s-D, -b-n-s, a+n-s+1]
        bb = [zero(T), -b-s, 1-D]
        return PaperGFormula("A", c*gamma(a+n+1), aa, bb, 2, 1, z)
    elseif row == :B
        aa = [1-s-D, -n-s-alpha]
        bb = [zero(T), -s-alpha, 1-D]
        return PaperGFormula("B", c, aa, bb, 1, 2, z)
    elseif row == :C
        aa = [1-s-D, -b-n-s, a+n-s+1, -s]
        bb = [zero(T), 1-D, -b-s, -s]
        return PaperGFormula("C", c*gamma(a+n+1), aa, bb, 1, 3, z)
    elseif row == :D
        aa = [1-s-D, a-s+1, a+b-s+1, -s]
        bb = [zero(T), -n-s, a+b+n-s+1, 1-D]
        return PaperGFormula("D", c*gamma(a+n+1), aa, bb, 1, 3, z)
    end
    throw(ArgumentError("row must be :A, :B, :C, or :D"))
end

end
