module Convexhull
import LinearAlgebra
import JuMP

import MathOptInterface
import NLopt


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket



import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots

function test()
    points = N -> [randn(2) for i in 1:N]
    v = points(30)
    # println(typeof(v))
    hull = convex_hull(v)
    # println(hull)
    typeof(hull), length(v), length(hull)

    p = plot([Singleton(vi) for vi in v])
    #println(p)
    plot!(p, VPolygon(hull), alpha=0.2)
end



function get_convexhull(socket::mySocket)
    #find a closed hull for the socket
    n = size(socket.points, 1)
    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n
        push!(ps, [socket.points[i].x, socket.points[i].y])
    end
    hull_ps = convex_hull(ps)
    hull = get_closed_socket(hull_ps)
    return hull
end

function get_closed_socket(ps)
    n = length(ps)
    socket = myType.mySocket()
	socket.points = Array{myType.Indexed_point}(undef, n)
	for i = 1:n

		socket.points[i] = myType.Indexed_point(i,
			ps[i][1], ps[i][2])
	end

	socket.edges = Array{myType.reference_edge}(undef, n)
	for i = 1:n
        next = (i+1)%n
        if next ==0
            next =n
        end
		socket.edges[i] = myType.reference_edge(i, i, next);

	end

    return socket

end

function get_index_set(index,n1)
    # get the point index for the pocket
    index_set = Array{Array{Int}}(undef,0)
    n = length(index)

    for i = 1:n

        next = (i+1)%n
        if next == 0
            next = n
        end

        if index[next] - index[i] == 1 || index[next] - index[i] == 1-n1
            continue
        end

        push!(index_set, [index[i], index[next]])
    end

    return index_set



end

function get_convexpocket(socket::mySocket)
    # return the pockets and hull

    n1 = length(socket.points)

    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n1
        push!(ps, [socket.points[i].x, socket.points[i].y])
    end

    hull_ps = convex_hull(ps)



    n2 = length(hull_ps)


    if n1 == n2
        error("no convexpocket")
    end




    index_hull = Array{Int}(undef, 0)
    index_pocket = Array{Int}(undef, 0)

    for i = 1:n1
        for p in hull_ps
            if ps[i][1]==p[1]&&ps[i][2]==p[2]
                push!(index_hull, i)
            else
                push!(index_pocket, i)
            end
        end
    end

    index_set = get_index_set(index_hull, n1)

    pockets = Array{mySocket}(undef, 0)

    #println(index_set)

    for id in index_set
        #println(id)
        if id[1]>id[2]
            error("no such socket")
        end

        p_pockets = Array{Indexed_point}(undef, 0)
        for i = id[1]:id[2]
            push!(p_pockets, socket.points[i])
        end
        push!(pockets, mySocket(p_pockets))
    end

    p_hull = Array{Indexed_point}(undef,0)

    for j in index_hull
        push!(p_hull, socket.points[j])
    end

    hull = mySocket(p_hull)




    return hull, pockets

end


function test_convexhull(i)
    peg,socket = Construct_Peg_Socket.read_from_file(i)
    hull = get_convexhull(socket)
    n = size(socket.points, 1)
    ps = Array{Array{Float64,1}}(undef, 0)
    for i = 1:n
        push!(ps, [socket.points[i].x, socket.points[i].y])
    end

    # println(ps)
    # println(hull)
    p = plot([Singleton(vi) for vi in ps])
    plot!(p, VPolygon(hull), alpha=0.2)
    #plot!(VPolygon(hull), alpha=0.2)

end

function test_pocket(i)
    peg,socket = Construct_Peg_Socket.read_from_file(i)
    hull, pockets = get_convexpocket(socket)





    for i = 1:length(pockets)
        c=(1, 1, 0)
        lw = 2.0
        Visualization.draw_sth(pockets[i],c, lw)
    end

    Visualization.draw_sth(hull)

end




#
end
