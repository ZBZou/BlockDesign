module Partial_order
import LinearAlgebra
import JuMP

import MathOptInterface
import NLopt


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket
import ..Insertion
import ..Convexhull

import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots

function get_poset(peg::myType.myPeg, socket::myType.mySocket, closed = false, convex =false)


    all_peg, c_p_s, model_list = Insertion.all_combinations(peg, socket, closed, convex)

    ms = []

    for i = 1:length(c_p_s)
        push!(ms, c_p_s[i])
    end



    for i = 1:length(model_list)

        for j = 1:length(model_list[i])
            push!(ms, model_list[i][j])
        end


    end



    if ms == []

        return nothing
    end

    n = length(ms)
    mo = length(ms[n])

    p_order = Array{Int}(undef, 0)






    for  i = 1:n

        if length(ms[i]) == mo
            push!(p_order, i)

        end

    end


    mo = mo -1

    if mo < 1
        return p_order
    else
        for i = mo:-1:1
            for j = 1:n

                if  length(ms[j]) == i

                    is_p = true

                    for k in p_order
                        if length(setdiff(ms[j], ms[k])) == 0
                            is_p = false
                            break
                        end
                    end

                    if is_p
                        push!(p_order, j)
                    end
                end
            end
        end
    end


    sort!(p_order)
    return p_order
end

function test_pd(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    closed = false
    convex = Insertion.is_convex(socket)
    p_order = get_poset(peg, socket, closed, convex)


end





end
