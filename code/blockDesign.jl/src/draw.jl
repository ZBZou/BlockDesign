module Draw

import Base.Iterators

using PyCall
using GraphPlot
using LightGraphs
using Compose
using Gadfly


import PyPlot
import JSON
import LightGraphs: vertices, edges, src, dst
import MetaGraphs: props
import GeometryTypes
import GeometryTypes: Point2f0
import CoordinateTransformations
import LinearAlgebra


import ..myType

import ..Construct_Peg_Socket

function graph()
    g = graphfamous("karate")
    gplot(g)
    #draw(PNG("mygraph.png", 8cm, 8cm), gplot(g))
end

function test()
    for i = 1:3
        x  = 1:21
        y1 = randn(21)
        y2 = randn(21)

        PyPlot.figure()
        PyPlot.axis("off")
        PyPlot.plot(x,y1,color="red",marker="*",linestyle="none")



        #PyPlot.show()
        #PyPlot.ion()

        PyPlot.figure()
        PyPlot.axis("off")
        PyPlot.plot(x,y2,color="blue",marker="*",linestyle="dashed")
        filename = "./anime/" * string(i) * ".png"
        PyPlot.savefig(filename)
    end
    #PyPlot.show()


end

function draw_sth(sth, c=(1, 0, 0), lw=2.0, title = "")


    p = Array{myType.Indexed_point}(undef, 0)

    for i = 1:length(sth.points)
        push!(p, sth.points[i])
    end

    push!(p, sth.points[1])



    #draw_polygon(p,c,lw)
    draw_line_seq(p, c, lw, title)

    #PyPlot.show()

end

function draw_open_sth(sth, c=(1, 0, 0), lw=2.0, title = "")
    p = Array{myType.Indexed_point}(undef, 0)

    for i = 1:length(sth.points)
        push!(p, sth.points[i])
    end





    #draw_polygon(p,c,lw)
    draw_line_seq(p, c, lw, title)

    #PyPlot.show()

end

function draw_anime(peg, socket, ori = [], m = [])

    if ori !=[]
        ori_peg, ori_socket = Construct_Peg_Socket.read_from_file(ori)
    end

    PyPlot.clf()
    PyPlot.axis("off")
    if ori !=[]
        Draw.draw_sth(ori_peg, (0,0,1), 3.)
        Draw.draw_open_sth(ori_socket, (0,0,0), 3.)
    end
    Draw.draw_sth(peg, (1,0,0), 3.)
    Draw.draw_open_sth(socket, (0.5,.4,.2), 3.)

    if m != []
        filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
        PyPlot.savefig(filename)
    end


end

function draw_peg_and_socket(peg, socket)
    Draw.draw_sth(peg, (0,0,1), 3.)
    Draw.draw_open_sth(socket, (0,0,0), 3.)
end

function draw_line_seq(points::Array{myType.Indexed_point}, c=(1, 0, 0), lw=2.0, title = "")
	# a sequence of points, not necessarily a closed shape
	len = size(points, 1)
	# the number of points


	for i = 1:len - 1
		x1 = points[i].x
		y1 = points[i].y

		x2 = points[i+1].x
		y2 = points[i+1].y



		PyPlot.plot([x1, x2], [y1, y2], color=c, linewidth=lw)
        PyPlot.axis("equal")
        PyPlot.axis("off")
        PyPlot.title(title)

	end



end

function draw_edges(egs::myType.mySocketEdges, c=(1, 0, 0), lw=2.0)
	# a sequence of points, not necessarily a closed shape
	n = length(egs.edges)
	# the number of points


	for i = 1:n - 1

        x1 = egs.edges[i].s.x
		y1 = egs.edges[i].s.y

		x2 = egs.edges[i].e.x
		y2 = egs.edges[i].e.y



		PyPlot.plot([x1, x2], [y1, y2], color=c, linewidth=lw)


	end



end

function draw_mlti_sth(socket::Array{myType.mySocket},peg, c=(1, 0, 0), lw=2.0)
    #PyPlot.figure()
    for i = 1:length(socket)
        p = Array{myType.Indexed_point}(undef, 0)

        for j = 1:length(socket[i].points)
            push!(p, socket[i].points[j])
        end

        draw_line_seq(p, (1,1,0), lw)
    end

    p = Array{myType.Indexed_point}(undef, 0)

    for i = 1:length(peg.points)
        push!(p, peg.points[i])
    end

    push!(p, peg.points[1])



    #draw_polygon(p,c,lw)
    draw_line_seq(p, c, lw)

    #PyPlot.show()
end

function draw_all_peg(all_peg::Array{myType.myPeg}, peg::myType.myPeg,
    socket::myType.mySocket)
    #PyPlot.figure()

    for i = 1: length(all_peg)
        PyPlot.figure()
        draw_open_sth(socket, (0,0,1), 2.0)
        draw_sth(peg, (1,0,1), 2.0)
        draw_sth(all_peg[i], (0,1,1), 2.0)
        # PyPlot.pause(5)
        # PyPlot.close()

    end


end


function draw_all_peg(all_peg::Array,  socket::myType.mySocket, op=true)
    #PyPlot.figure()
    if op == false
        for i = 1: length(all_peg)
            PyPlot.figure()
            PyPlot.axis("off")
            draw_sth(socket, (0,0,0), 2.0)
            draw_sth(all_peg[i], (0,0,1), 2.0)
        # PyPlot.pause(5)
        # PyPlot.close()

        end
    else
        for i = 1: length(all_peg)
            PyPlot.figure()
            PyPlot.axis("off")
            draw_open_sth(socket, (0,0,0), 2.0)
            draw_sth(all_peg[i], (0,0,1), 2.0)
        # PyPlot.pause(5)
        # PyPlot.close()

        end
    end


end

function draw_all_peg_with_horizontal_line(all_peg::Array,  socket::myType.mySocket, title_list = "", op=true)
    #PyPlot.figure()
    n = length(socket.points)
    if op == false
        for i = 1: length(all_peg)
            if title_list!=""
                title = title_list[i]
            else
                title = ""
            end
            PyPlot.figure()
            PyPlot.axis("off")
            draw_sth(socket, (0,0,0), 1.0, title)
            draw_sth(all_peg[i], (0,0,1), 1.0, title)
            p1 = myType.Indexed_point(0, -4., socket.points[1].y)
            p2 = myType.Indexed_point(n+1, 4., socket.points[n].y)
            draw_line_seq([p1, socket.points[1]], (0,0,0), 1.0, title)
            draw_line_seq([p2, socket.points[n]], (0,0,0), 1.0, title)

        PyPlot.pause(2)
        PyPlot.close()

        end
    else
        for i = 1: length(all_peg)
            if title_list!=""
                title = title_list[i]
            else
                title = ""
            end
            PyPlot.figure()
            PyPlot.axis("off")
            draw_open_sth(socket, (0,0,0), 1.0, title)
            draw_sth(all_peg[i], (0,0,1), 1.0, title)
            p1 = myType.Indexed_point(0, -4., socket.points[1].y)
            p2 = myType.Indexed_point(n+1, 4., socket.points[n].y)
            draw_line_seq([p1, socket.points[1]], (0,0,0), 1.0, title)
            draw_line_seq([p2, socket.points[n]], (0,0,0), 1.0, title)
        PyPlot.pause(2)
        PyPlot.close()

        end
    end


end

function test2(i)
    peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
    draw_sth(socket)
    draw_sth(peg, (0,0,1))
end

function draw_peg_socket(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    draw_sth(socket)
    draw_sth(peg, (0,0,1))
    
    # PyPlot.plot([x1, x2], [y1, y2])
    
end

function norm(v)
    sum = 0
    for i = 1:length(v)
        sum += v[i]*v[i]
    end
    return sqrt(sum)
end

function dist(p1, p2)
    dist = sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2)
    return dist
end

function test_dist(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    ps = socket.points
    # print(length(ps))
    # 
    l1 = [ps[1], ps[2]]
    l2 = [ps[2], ps[3]]
    l3 = [ps[3], ps[4]]
    l4 = [ps[4], ps[5]]
    l5 = [ps[5], ps[6]]
    l6 = [ps[6], ps[7]]
    
    p1 = peg.points[1]
    p2 = peg.points[2]
    p3 = peg.points[3]
    p4 = peg.points[4]
    p5 = peg.points[5]
    p6 = peg.points[6]
    # println("ps: ", p1, ", ", p2, ", ", p3, ", ", p4, ", ", p5)
    
    
    v1 = [ps[2].x - ps[1].x, ps[2].y - ps[1].y]
    v2 = [ps[3].x - ps[2].x, ps[3].y - ps[2].y]
    v3 = [ps[4].x - ps[3].x, ps[4].y - ps[3].y]
    v4 = [ps[5].x - ps[4].x, ps[5].y - ps[4].y]
    v5 = [ps[6].x - ps[5].x, ps[6].y - ps[5].y]
    v6 = [ps[7].x - ps[6].x, ps[7].y - ps[6].y]
    
    v1 = v1 / norm(v1)
    v2 = v2 / norm(v2)
    v3 = v3 / norm(v3)
    v4 = v4 / norm(v4)
    v5 = v5 / norm(v5)
    v6 = v6 / norm(v6)
    
    n1 = [-v1[2], v1[1]]
    n2 = [-v2[2], v2[1]]
    n3 = [-v3[2], v3[1]]
    n4 = [-v4[2], v4[1]]
    n5 = [-v5[2], v5[1]]
    n6 = [-v6[2], v6[1]]
    
    lp1 = [p1.x - ps[1].x, p1.y - ps[1].y]
    lp2 = [p2.x - ps[2].x, p2.y - ps[2].y]
    lp3 = [p3.x - ps[3].x, p3.y - ps[3].y]
    lp4 = [p4.x - ps[4].x, p4.y - ps[4].y]
    lp5 = [p5.x - ps[5].x, p5.y - ps[5].y]
    lp6 = [p6.x - ps[6].x, p6.y - ps[6].y]
    
    # lp4 = [p3.x - ps[4].x, p3.y - ps[4].y]
    # lp5 = [p4.x - ps[5].x, p4.y - ps[5].y]
    # lp6 = [p5.x - ps[6].x, p5.y - ps[6].y]
    
    temp = v1[1]*lp1[1] + v1[2]*lp1[2]
    np1 = [ps[1].x + temp * v1[1], ps[1].y + temp * v1[2]]
    
    temp = v2[1]*lp2[1] + v2[2]*lp2[2]
    np2 = [ps[2].x + temp * v2[1], ps[2].y + temp * v2[2]]
    
    temp = v3[1]*lp3[1] + v3[2]*lp3[2]
    np3 = [ps[3].x + temp * v3[1], ps[3].y + temp * v3[2]]
    
    temp = v4[1]*lp4[1] + v4[2]*lp4[2]
    np4 = [ps[4].x + temp * v4[1], ps[4].y + temp * v4[2]]
    
    temp = v5[1]*lp5[1] + v5[2]*lp5[2]
    np5 = [ps[5].x + temp * v5[1], ps[5].y + temp * v5[2]]
    
    temp = v6[1]*lp6[1] + v6[2]*lp6[2]
    np6 = [ps[6].x + temp * v6[1], ps[6].y + temp * v6[2]]
    
    # currently, 0.1 error
    # produce different error 
    # 
    # 
    # 0.2 
    # 
    scale = 0.2
    d02p1 = [np1[1]+scale*n1[1], np1[2]+scale*n1[2]]
    d02p2 = [np2[1]+scale*n2[1], np2[2]+scale*n2[2]]
    d02p3 = [np3[1]+scale*n3[1], np3[2]+scale*n3[2]]
    d02p4 = [np4[1]+scale*n4[1], np4[2]+scale*n4[2]]
    d02p5 = [np5[1]+scale*n5[1], np5[2]+scale*n5[2]]
    d02p6 = [np6[1]+scale*n6[1], np6[2]+scale*n6[2]]
    d02p6 = []
    println("dist 0.2: ", d02p1, ", ", d02p2, ", ", d02p3, ", ", d02p4, ", ", d02p5, ", ", d02p6)
    
    
    # 
    # 0.15
    # 
    scale = 0.15
    d015p1 = [np1[1]+scale*n1[1], np1[2]+scale*n1[2]]
    d015p2 = [np2[1]+scale*n2[1], np2[2]+scale*n2[2]]
    d015p3 = [np3[1]+scale*n3[1], np3[2]+scale*n3[2]]
    d015p4 = [np4[1]+scale*n4[1], np4[2]+scale*n4[2]]
    d015p5 = [np5[1]+scale*n5[1], np5[2]+scale*n5[2]]
    d015p6 = [np6[1]+scale*n6[1], np6[2]+scale*n6[2]]
    d015p6 = []
    println("dist 0.15: ", d015p1, ", ", d015p2, ", ", d015p3, ", ", d015p4, ", ", d015p5, ", ", d015p6)
    
    
    # 
    # 
    # 0.05 
    # 
    # 
    scale = 0.05
    d005p1 = [np1[1]+scale*n1[1], np1[2]+scale*n1[2]]
    d005p2 = [np2[1]+scale*n2[1], np2[2]+scale*n2[2]]
    d005p3 = [np3[1]+scale*n3[1], np3[2]+scale*n3[2]]
    d005p4 = [np4[1]+scale*n4[1], np4[2]+scale*n4[2]]
    d005p5 = [np5[1]+scale*n5[1], np5[2]+scale*n5[2]]
    d005p6 = [np6[1]+scale*n6[1], np6[2]+scale*n6[2]]
    d005p6 = []
    println("dist 0.05: ", d005p1, ", ", d005p2, ", ", d005p3, ", ", d005p4, ", ", d005p5, ", ", d005p6)
    # 
    # 
    # 0.03
    # 
    # scale = 0.03
    # d003p1 = [np1[1]+scale*n1[1], np1[2]+scale*n1[2]]
    # d003p2 = [np2[1]+scale*n2[1], np2[2]+scale*n2[2]]
    # d003p3 = [np3[1]+scale*n3[1], np3[2]+scale*n3[2]]
    # d003p4 = [np4[1]+scale*n4[1], np4[2]+scale*n4[2]]
    # d003p5 = [np5[1]+scale*n5[1], np5[2]+scale*n5[2]]
    # # d003p6 = [np6[1]+scale*n6[1], np6[2]+scale*n6[2]]
    # 
    # print("dist 0.03: ", d003p1, ", ", d003p2, ", ", d003p3, ", ", d003p4, ", ", d003p5, ", ", d003p6)
    # 
    # 
    # 
    # 0.01
    # 
    scale = 0.005
    d001p1 = [np1[1]+scale*n1[1], np1[2]+scale*n1[2]]
    d001p2 = [np2[1]+scale*n2[1], np2[2]+scale*n2[2]]
    d001p3 = [np3[1]+scale*n3[1], np3[2]+scale*n3[2]]
    d001p4 = [np4[1]+scale*n4[1], np4[2]+scale*n4[2]]
    d001p5 = [np5[1]+scale*n5[1], np5[2]+scale*n5[2]]
    d001p6 = [np6[1]+scale*n6[1], np6[2]+scale*n6[2]]
    d001p6 = []
    println("dist 0.01: ", d001p1, ", ", d001p2, ", ", d001p3, ", ", d001p4, ", ", d001p5, ", ", d001p6)
    
    
    
end



end
