module Main

import LinearAlgebra
import JuMP
import Luxor

import MathOptInterface
import NLopt


import.. CONST
import.. ConcaveTest
import.. Construct_Peg_Socket
import.. Draw
import.. Kinematics
import.. Constraints
import.. ContactMode
import.. CMT
import.. Convexhull
import.. Adjust
import.. Insertion
import.. Partial_order
import.. Rotation
import.. CMT_graph
import.. Insertion_gradient
import.. Stability_gradient
import.. Opjoint


import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots

# There are three main steps: 1) eliminate_undesired_sink; 2) peg_optimizing and 3) socket_optimize.
# dth and dt both are optimization step size.
function optimizing_joint(joint_id, cmt_and_sinks_figures = false, anime = false, rec_error = false, rec_stab = false)
    peg, socket = Construct_Peg_Socket.read_from_file(joint_id)
    c_ids = Construct_Peg_Socket.get_contact_ids(joint_id)
    peg = Adjust.get_adj_peg(peg, socket, 0.)
    dt1 = [0.2, 0.1]
    dt2 = [0.1, 0.1]
    dt3 = [0.1, 0.1]
    dt4 = [0.1, 0.1]
    dt5 = [0.1, 0.1]
    dt6 = [0.1, 0.1]
    dt7 = [0.1, 0.1]
    dt8 = [0.1, 0.1]
    dt9 = [0.1, 0.1]
    dt10 = [0.1, 0.1]
    dt11 = [0.1, 0.1]
    dt12 = [0.1, 0.1]
    # dth = [dt1, dt2, dt3, dt4]
    # dth = [dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8]
    dth = [dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8, dt9, dt10, dt11, dt12]
    
    dt = 0.01
    pic_count = 0
    before_error = 0
    after_error = 0
    before_stab = 0
    after_stab = 0

    # if cmt_and_sinks_figures
    #     CMT_graph.cmt_graph(peg, socket, c_ids)
    # end

    println("joint_id: ", joint_id)

    new_peg1, new_socket1, pic_count = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids, joint_id, pic_count, dt)
    new_peg2, new_socket2, pic_count = Insertion_gradient.peg_optimizing(new_peg1, new_socket1, c_ids, joint_id, pic_count, dth[joint_id][1])
    new_peg3, new_socket3, pic_count = Stability_gradient.socket_optimize(new_peg2, new_socket2, c_ids, joint_id, pic_count, dth[joint_id][2])
    

    # println(new_peg3, new_socket3)
    if cmt_and_sinks_figures
        CMT_graph.cmt_graph(new_peg3, new_socket3, c_ids)
    end
    
    # rec_error = true 
    rec_stab = true

    if rec_error

        sinks = CMT_graph.get_sinks(new_peg1,new_socket1,c_ids)
        desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
        before_error = CMT_graph.sink_with_maximum_error_output_error(new_peg1,new_socket1,desired_sinks)
        new_sinks = CMT_graph.get_sinks(new_peg3,new_socket3,c_ids)
        new_desired_sinks = CMT_graph.desired_sinks(c_ids,new_sinks)
        after_error = CMT_graph.sink_with_maximum_error_output_error(new_peg3,new_socket3,new_desired_sinks)
        println("to print insertion error")
        fname = "insertion_error.txt"
        open(fname,"a") do io
             println(io,"index =",joint_id)
             println(io,"before_error =",before_error)
             println(io,"after_error =",after_error)
        end

    end

    if rec_stab
        sinks = CMT_graph.get_sinks(new_peg1,new_socket1,c_ids)
        desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
        before_stab = Stability_gradient.total_stab(new_peg1, new_socket1, desired_sinks)
        # before_stab, before_sink = Stability_gradient.sink_with_least_stab(new_peg1, new_socket1, desired_sinks)

        new_sinks = CMT_graph.get_sinks(new_peg3,new_socket3,c_ids)
        new_desired_sinks = CMT_graph.desired_sinks(c_ids,new_sinks)
        after_stab = Stability_gradient.total_stab(new_peg3, new_socket3, new_desired_sinks)

        # after_stab, after_sink = Stability_gradient.sink_with_least_stab(new_peg3, new_socket3, new_desired_sinks)
        println("to print stability")
        fname = "stability.txt"
        open(fname,"a") do io
             println(io,"index =",joint_id)
             println(io,"before_stab =",before_stab)
             println(io,"after_stab =",after_stab)
        end
    end

    jointfile = "joint.txt"

    for  i = 1: length(new_peg1.points)
        open(jointfile,"a") do io

             println(io,"peg:[index, x, y] =", [i, new_peg1.points[i].x, new_peg1.points[i].y])

        end
    end

    for  j = 1: length(new_socket1.points)
        open(jointfile,"a") do io

             println(io,"socket:[index, x, y] =", [j, new_socket1.points[j].x, new_socket1.points[j].y])

        end
    end

    for  i = 1: length(new_peg3.points)
        open(jointfile,"a") do io

             println(io,"peg:[index, x, y] =", [i, new_peg3.points[i].x, new_peg3.points[i].y])

        end
    end

    for  j = 1: length(new_socket3.points)
        open(jointfile,"a") do io

             println(io,"socket:[index, x, y] =", [j, new_socket3.points[j].x, new_socket3.points[j].y])

        end
    end


end

function test1(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
    t = Adjust.compute_t(peg, socket)
    new_socket = Rotation.adj_upsocket(socket, 0.5)
    new_peg = Adjust.get_peg_from_t(peg, new_socket, t)
    c_p_s, cm_ids, c_lists = ContactMode.find_valid_three_pairs_cms(new_peg, new_socket, c_ids)
    sinks = CMT_graph.get_sinks(new_peg, new_socket, c_ids)
    sink = CMT_graph.sink_with_maximum_error(new_peg, new_socket, sinks)
    println("cm_ids: ", cm_ids)
    println("sinks: ", sink)
    Draw.draw_sth(new_peg, (0,0,1))
    Draw.draw_open_sth(new_socket)
end


function test2(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
    Insertion_gradient.peg_optimizing(peg, socket, c_ids)
end

function test3(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
    peg = Adjust.get_adj_peg(peg, socket, [0.5, 0.5])
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
    desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
    # println(sinks)
    sink = CMT_graph.sink_with_maximum_error(peg,socket,desired_sinks)
	k = Insertion_gradient.peg_gradient_to_minimize_insertion_error(peg, socket, sink)
end

function test4(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    # peg = Adjust.get_adj_peg(peg, socket, [0.5, 0.5])
    t = [0.41011235955056174, 0.9459405940594061]

    peg = Adjust.get_peg_from_t(peg, socket, t)
    peg = Adjust.get_adj_peg(peg, socket, 0.)
    cps, pegs = ContactMode.find_all_pairs(peg, socket, c_ids)
    println(cps)
    cm = myType.cm_id([cps[3], cps[6], cps[11]])
    # println(cm)
    result, c_list= ContactMode.is_cm_valid(peg, socket, cm)

end

function test5(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
    # CMT_graph.cmt_graph(peg, socket, c_ids)
	new_peg, new_socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids, 0.01)
	new_peg, new_socket = Insertion_gradient.peg_optimizing(new_peg, new_socket, c_ids)
    # CMT_graph.cmt_graph(peg, socket, c_ids)
    t = Adjust.compute_t(new_peg, new_socket)
    v1 = [socket.points[2].x - socket.points[1].x, socket.points[2].y - socket.points[1].y, 0]
    v2 = [new_socket.points[2].x - new_socket.points[1].x, new_socket.points[2].y - new_socket.points[1].y, 0]
    ag = Constraints.compute_signed_angle(v1, v2)
    println(t, ag)
end

function test6(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    sinks = CMT_graph.get_sinks(peg,socket,c_ids)
    desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
    # println(sinks)
    # sink = CMT_graph.sink_with_maximum_error(peg,socket,desired_sinks)
    println(desired_sinks)
    stab_1 = Stability_gradient.sink_stability(peg, socket, desired_sinks[1])
    stab_2 = Stability_gradient.sink_stability(peg, socket, desired_sinks[2])
    stab_3 = Stability_gradient.sink_stability(peg, socket, desired_sinks[3])
    stab_4 = Stability_gradient.sink_stability(peg, socket, desired_sinks[4])
    println([stab_1, stab_2, stab_3, stab_4])
	# peg, socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids, 0.01)
    peg, socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids)
    sinks = CMT_graph.get_sinks(peg,socket,c_ids)
    desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
    # println(sinks)
    # sink = CMT_graph.sink_with_maximum_error(peg,socket,desired_sinks)
    println(desired_sinks)
    stab_1 = Stability_gradient.sink_stability(peg, socket, desired_sinks[1])
    stab_2 = Stability_gradient.sink_stability(peg, socket, desired_sinks[2])
    stab_3 = Stability_gradient.sink_stability(peg, socket, desired_sinks[3])
    stab_4 = Stability_gradient.sink_stability(peg, socket, desired_sinks[4])
    println([stab_1, stab_2, stab_3, stab_4])

    # cps, pegs = ContactMode.find_all_pairs(peg, socket, c_ids)
    # println(cps)
    # cm = myType.cm_id([cps[3], cps[6], cps[11]])
    # # println(cm)
    # result, c_list= ContactMode.is_cm_valid(peg, socket, cm)

end

function test7(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    t = [0.1501123595505622, 0.8559405940594067]
    ang = 0.33
    socket1 = Rotation.adj_upsocket(socket, ang)
    peg1 = Adjust.get_peg_from_t(peg, socket1, t)
    sinks1 = CMT_graph.get_sinks(peg1,socket1,c_ids)
    desired_sinks1 = CMT_graph.desired_sinks(c_ids,sinks1)
    stab1 = Stability_gradient.total_stab(peg1, socket1, desired_sinks1)

    dth = 0.1
    k1 = -1
    r = [0,0]
    index = 1
    socket2 = Rotation.rotate_socket(socket1, r, index, dth*k1)
    peg2 = Adjust.get_adj_peg(peg1, socket2, 0.)

    stab2 = Stability_gradient.total_stab(peg2, socket2, desired_sinks1)

    # stab2 = Stability_gradient.sink_stability(peg2, socket2, sink1)

    println([stab1, stab2])


    # Draw.draw_sth(new_peg)
    Draw.draw_open_sth(socket1)
    Draw.draw_open_sth(socket2, (0,0,1))
    # sinks = CMT_graph.get_sinks(new_peg,new_socket,c_ids)
    # r, index, k = Stability_gradient.stab_gradient(new_peg, new_socket, c_ids)
end

function test8(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    t = [0.1501123595505622, 0.8559405940594067]
    ang = 0.33
    new_socket = Rotation.adj_upsocket(socket, ang)
    new_peg = Adjust.get_peg_from_t(peg, new_socket, t)

    # Draw.draw_sth(new_peg)
    # Draw.draw_open_sth(new_socket)
    r, index, k = Stability_gradient.stab_gradient(new_peg, new_socket, c_ids)
end

function test9(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    t = [0.1501123595505622, 0.8559405940594067]
    ang = 0.33
    new_socket = Rotation.adj_upsocket(socket, ang)
    new_peg = Adjust.get_peg_from_t(peg, new_socket, t)
    Stability_gradient.socket_optimize(new_peg, new_socket, c_ids)


end

function test10(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    println("c_ids: ", c_ids)
    PyPlot.figure()
    Draw.draw_sth(peg, (0,0,1), 3.)
    Draw.draw_open_sth(socket, (0,0,0), 3.)
    PyPlot.figure()
    CMT_graph.cmt_graph(peg, socket, c_ids)

    PyPlot.figure()
	# new_peg, new_socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids, 0.01)
    new_peg, new_socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids)
	new_peg, new_socket = Insertion_gradient.peg_optimizing(new_peg, new_socket, c_ids)
    Stability_gradient.socket_optimize(new_peg, new_socket, c_ids)
    PyPlot.figure()
    CMT_graph.cmt_graph(peg, socket, c_ids)


end


end  # module Main
