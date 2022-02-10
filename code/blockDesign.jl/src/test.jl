module Test

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
import.. Main
import.. Test

import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots

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
    # peg = Adjust.get_adj_peg(peg, socket, [0.5, 0.5])
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
	peg, socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids, 0.01)
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
    PyPlot.figure()
    Draw.draw_sth(peg, (0,0,1), 3.)
    Draw.draw_open_sth(socket, (0,0,0), 3.)
    PyPlot.figure()
    CMT_graph.cmt_graph(peg, socket, c_ids)

    PyPlot.figure()
	new_peg, new_socket = Insertion_gradient.eliminate_undesired_sink(peg, socket, c_ids, 0.01)
	new_peg, new_socket = Insertion_gradient.peg_optimizing(new_peg, new_socket, c_ids)
    Stability_gradient.socket_optimize(new_peg, new_socket, c_ids)
    PyPlot.figure()
    CMT_graph.cmt_graph(peg, socket, c_ids)


end

function test11(i)
    cmt_and_sinks_figures = true
    anime = true
    rec_error = true
    rec_stab = true

    Main.optimizing_joint(i,cmt_and_sinks_figures, anime, rec_error, rec_stab)
end

function test12(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    println(c_ids)
    CMT_graph.cmt_graph(peg, socket, c_ids)
end

function test13(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    new_peg = Adjust.get_adj_peg(peg, socket, 0.)

    # Draw.draw_peg_and_socket(peg, socket)
    # Draw.draw_sth(new_peg)
    # println(new_peg)

end

function test14(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    new_peg = Adjust.get_adj_peg(peg, socket, 0.)
    # Draw.draw_peg_and_socket(peg, socket)
    # Draw.draw_sth(new_peg)
    es = myType.mySocketEdges(socket)
    cp = myType.contact_pair(1, peg.points[1], es.edges[1])
    result, config_list = ContactMode.is_contact_pair_valid(peg, socket, cp)
    # c_p_s, cm_ids, c_lists = ContactMode.find_valid_three_pairs_cms(peg, socket, c_ids)
    # cm = myType.cm_id(c_p_s[1], c_p_s, c_p_s])

    # cp1 = myType.contact_id(1, 1)
    # cp2 = myType.contact_id(3, 4)
    # cp3 = myType.contact_id(4, 5)
    # cm_id = myType.cm_id([cp1, cp2, cp3])
    #
    # result, c_list= ContactMode.is_cm_valid(peg, socket, cm_id)
    # pegs = Kinematics.get_pegs_from_configs(peg,c_lists)
    #
    # Draw.draw_all_peg(pegs, socket)
end

end  # module Test
