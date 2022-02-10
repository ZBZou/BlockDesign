module ConcaveTest
import ..myType
using LazySets
using Plots

function test()
    points = N -> [randn(2) for i in 1:N]
    v = points(30)
    #println(typeof(v))
    hull = convex_hull(v)
    # println(hull)
    # typeof(hull), length(v), length(hull)

    p = plot([Singleton(vi) for vi in v])
    #println(p)
    plot!(p, VPolygon(hull), alpha=0.2)
end

function is_concave(socket::myType.mySocket)
    n = length(socket.points)
    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n
        push!(ps, [socket.points[i].x, socket.points[i].y])
    end
    hull_ps = convex_hull(ps)
    if length(hull_ps) == n
        return false
    else
        return true
    end

end

function is_concave(peg::myType.myPeg)
    n = length(peg.points)
    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n
        push!(ps, [peg.points[i].x, peg.points[i].y])
    end
    hull_ps = convex_hull(ps)
    if length(hull_ps) == n
        return false
    else
        return true
    end

end

function concave_points(peg::myType.myPeg)
    id = []
    n = length(peg.points)
    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n
        push!(ps, [peg.points[i].x, peg.points[i].y])
    end
    hull_ps = convex_hull(ps)

    for i = 1:n
        if !(ps[i] in hull_ps)
             push!(id, i)
         end
     end

     return id
end

function concave_points(socket::myType.mySocket)
    id = []
    n = length(socket.points)
    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n
        push!(ps, [socket.points[i].x, socket.points[i].y])
    end
    hull_ps = convex_hull(ps)

    for i = 1:n
        if !(ps[i] in hull_ps)
             push!(id, i)
         end
     end

     return id
end

end  # module Concave_test
