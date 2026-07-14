const POINTS_PER_FORMULA = 100
const OSCILLATORY_OPERATOR_POINTS = 20
const SEED = 0x3872026

relative_error(x, y) = abs(x-y) / max(abs(x), abs(y), 1e-12)
distance_to_integer(x) = abs(x - round(x))

function safely_nonconfluent(f)
    left_b = f.b[1:f.m]
    left_a = f.a[1:f.n]
    lower_sep = all(distance_to_integer(left_b[i]-left_b[j]) > 0.075
                    for i in eachindex(left_b) for j in 1:(i-1))
    upper_sep = all(distance_to_integer(left_a[i]-left_a[j]) > 0.075
                    for i in eachindex(left_a) for j in 1:(i-1))
    lower_sep && (!(length(f.a) == length(f.b) && f.z > 1) || upper_sep)
end

function classical_laplacian_highd(row, r; d, ell, n, a, b, alpha,
                                   step=2e-4)
    radial(x) = higher_source_radial(row,x;n,a,b,alpha)
    first = (radial(r-2step)-8radial(r-step)+8radial(r+step)-radial(r+2step)) /
            (12step)
    second = (-radial(r+2step)+16radial(r+step)-30radial(r)+
              16radial(r-step)-radial(r-2step)) / (12step^2)
    -r^ell*(second+(d+2ell-1)*first/r)
end
