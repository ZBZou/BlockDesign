module Insertion_gradient


using ..myType
using Compose

import LightGraphs
import GraphPlot
import GraphRecipes
import LinearAlgebra
import JuMP
import Cairo
import Fontconfig

import MathOptInterface
import NLopt


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket
import ..Convexhull
import ..Draw
import ..CMT
import ..Insertion
import ..Insertion
import ..CMT_graph
import ..Adjust
import ..ContactMode
import ..Rotation

import PyPlot
import Plots

function friction_cone(socket)
	es = myType.mySocketEdges(socket)
	n = length(es.edges)
	fn = []
	for i = 1:n
		fn_i = outward_normal(es.edges[i].s, es.edges[i].e)
		push!(fn, fn_i)
	end

	fcone = []
	for i = 1:n
		v1 = [fn[i][1], fn[i][2], 0]
		v2 = [0,-1,0]
		thi = Insertion.compute_signed_angle(v2, v1)
		fci = [thi - atan(CONST.friction_coeff), thi + atan(CONST.friction_coeff)]
		push!(fcone, fci)
	end
	return fcone
end

function force_analysis_for_sink(socket, peg, sink, theta)
	F = Insertion.th_get_F(theta, 2)
	Fs = []
	n = length(sink)
	es = myType.mySocketEdges(socket)
	for i = 1:n
        target_p = sink[i].p_id
        F_all = deepcopy(F)
        for j = 1:n
			if target_p == sink[j].p_id
				continue
			end

			F_j = deepcopy(F)
			cp_j = myType.contact_pair(0, peg.points[sink[j].p_id], es.edges[sink[j].e_id])
			F_pj, Fnj= Insertion.force_between_pair(cp_j, peg, socket, F_j)

			if Fnj>0
				F_all[1] += F_pj[1]
                F_all[2] += F_pj[2]
			end
		end
		v1 = [F_all[1], F_all[2], 0]
		v2 = [0,-1,0]
		thi = Insertion.compute_signed_angle(v2, v1)
		push!(Fs, thi)
    end
	return Fs
end

function outward_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(-pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function all_close(a, b)
	if abs(a - b) < 0.001
		return true
	end
	return false;

end

function peg_gradient_to_minimize_insertion_error(peg, socket, sink, delta_t = 0.1, limit_idex=[])
	cm = myType.cm_id(sink)
	th = abs(CMT_graph.insertion_error(peg, socket, cm))
	n_t = Int(floor(length(peg.points)/2))
	dt = []
	for i = 1:n_t
		push!(dt, delta_t)
	end

	g = zeros(n_t)
	k_new = zeros(n_t)
	k_id = Array{Array}(undef, 0)

	limit_t = Adjust.check_t_limit(peg, socket)
	# Draw.draw_open_sth(socket)
	# Draw.draw_sth(peg)
	for i = 1:n_t
		peg_i = deepcopy(peg)
		if limit_t[i] || i in limit_idex
			# println("limit_t: ", i)
			continue
		end

		g1 = deepcopy(g)
		g2 = deepcopy(g)
		g1[i] = 1
		g2[i] = -1

		peg1 = Adjust.get_adj_peg(peg_i, socket, g1.*dt)
		peg2 = Adjust.get_adj_peg(peg_i, socket, g2.*dt)



		th1= CMT_graph.insertion_error(peg1, socket, cm)
		th2= CMT_graph.insertion_error(peg2, socket, cm)
		# show([th1, th2])
		if th1 == nothing
			th1 = 100
		end
		if th2 == nothing
			th2 = 100
		end

		th1= abs(th1)
		th2= abs(th2)
		# println([i, th, th1, th2])

		if th <= th1 && th <= th2

			push!(k_id, [th, i, 0])
		elseif th1 < th2

			push!(k_id, [th1, i, 1])
		else

			push!(k_id, [th2, i, -1])
		end
	end

	if k_id != []
		sort!(k_id, by = x->x[1])
		id = Int(k_id[1][2])
		k_new[id] = k_id[1][3]
	end

	# println("peg k_id : ", k_id)
	# println("k_new to return: ", k_new)
	return k_new
end

function peg_optimizing(peg, socket, c_ids, i = [], pic_count = 0, dt = 0.1, k_dt=0.01)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
	sink = CMT_graph.sink_with_maximum_error(peg,socket,desired_sinks)
	if sink == nothing || length(sink) == 0
		return peg, socket, pic_count
	end
	k = peg_gradient_to_minimize_insertion_error(peg, socket, sink, dt)
	# println(k)
	output_pic_count = pic_count
	new_k = deepcopy(k)
	new_peg = deepcopy(peg)
	new_socket = deepcopy(socket)
	output_peg = deepcopy(new_peg)
	output_socket = deepcopy(new_socket)
	# println("peg_optimizing: ", dt, ", ", k, ", ", k_dt)
	while unique!(new_k) != [0]
		new_peg = Adjust.get_adj_peg(new_peg, new_socket, k*k_dt)
		sinks = CMT_graph.get_sinks(peg,socket,c_ids)
		# stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
		# if stuck_sink!= []
		# 	new_peg, new_socket, pic_count = eliminate_undesired_sink(peg, socket, c_ids, i, pic_count)
		# end
		k = peg_gradient_to_minimize_insertion_error(new_peg, new_socket, sink, dt)
		Draw.draw_anime(new_peg, new_socket, i, pic_count)
		new_k = deepcopy(k)
		# PyPlot.clf()
		# PyPlot.axis("off")
		# Draw.draw_sth(new_peg, (1,0,0), 3.)
		# Draw.draw_open_sth(new_socket, (0.5,.4,.2), 3.)
		output_peg = deepcopy(new_peg)
		output_socket = deepcopy(new_socket)
		output_pic_count = pic_count
		pic_count += 1
	end

	return output_peg, output_socket, output_pic_count
end

function test_peg_gradient(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
	k = peg_gradient_to_minimize_insertion_error(peg, socket, stuck_sink[1])
end

function edge_rotate(socket::myType.mySocket, th::Float64, e_id::Int, p_c::Array{Float64,1})

    n = length(socket.points)

    p_s = socket.points[socket.edges[e_id].s]
    p_e = socket.points[socket.edges[e_id].e]

    R = Kinematics.rotation_matrix_2D(th)

    ps = [p_s.x, p_s.y]
    pe = [p_e.x, p_e.y]

    next = (e_id+1)%n
    if next == 0
        next = n
    end

    pn = [socket.points[socket.edges[next].e].x, socket.points[socket.edges[next].e].y]

    pe_new = p_c + R*(pe-p_c)
    ps_new = p_c + R*(ps-p_c)



    l1 = get_line(ps_new,pe_new)
    l2 = get_line(pe,pn)


    pe_new = get_intersection(l1,l2)


    new_socket = deepcopy(socket)

    vn = pe_new - pe

    new_socket.points[socket.edges[e_id].s].x = ps_new[1]
    new_socket.points[socket.edges[e_id].s].y = ps_new[2]

    for i = e_id+1:n
        new_socket.points[i].x = socket.points[i].x + vn[1]
        new_socket.points[i].y = socket.points[i].y + vn[2]
    end

    return new_socket

end


function edge_rotate(egs::myType.mySocketEdges, th::Float64, e_id::Int, p_c::Array{Float64,1})

    n = length(egs.edges)

    p_s = egs.edges[e_id].s
    p_e = egs.edges[e_id].e

    R = Kinematics.rotation_matrix_2D(th)

    ps = [p_s.x, p_s.y]
    pe = [p_e.x, p_e.y]

    pe_new = p_c + R*(pe-p_c)
    ps_new = p_c + R*(ps-p_c)

    new_egs = deepcopy(egs)

    new_egs.edges[e_id].s.x = ps_new[1]
    new_egs.edges[e_id].s.y = ps_new[2]
    new_egs.edges[e_id].e.x = pe_new[1]
    new_egs.edges[e_id].e.y = pe_new[2]



    return new_egs

end


function get_line(p1::Array{Float64,1}, p2::Array{Float64,1})
    A = (p1[2] - p2[2])
    B = (p2[1] - p1[1])
    C = (p1[1]*p2[2] - p2[1]*p1[2])
    return [A, B, -C]
end

function get_line(p1::Indexed_point, p2::Indexed_point)
    A = (p1.y - p2.y)
    B = (p2.x - p1.x)
    C = (p1.x*p2.y - p2.x*p1.y)
    return [A, B, -C]
end

function get_intersection(L1::Array{Float64,1}, L2::Array{Float64,1})
    D  = L1[1] * L2[2] - L1[2] * L2[1]
    Dx = L1[3] * L2[2] - L1[2] * L2[3]
    Dy = L1[1] * L2[3] - L1[3] * L2[1]
    if D != 0
        x = Dx / D
        y = Dy / D
        return [x,y]
    else
        error("no intersection")
    end
end

function cp_not_move(c_p::myType.contact_pair, peg::myType.myPeg,
				socket::myType.mySocket, config::myType.configSE2,
				F::Array, mode::Int)

	# rotate edge untill the contact points of CP not moving
    new_peg = Kinematics.get_peg_from_config(peg, config)
    if mode == 0
        px = new_peg.points[c_p.c_p.index].x
        py = new_peg.points[c_p.c_p.index].y
    elseif mode == 1
        px = socket.points[c_p.c_e.index].x
        py = socket.points[c_p.c_e.index].y
    end

    egs = myType.mySocketEdges(socket)

    for th = -pi:pi/90:pi

        new_egs = edge_rotate(egs, th, c_p.c_e.index, [px, py])
        #Draw.draw_open_sth(new_socket)
        F_p, Fn = Insertion.force_between_pair(c_p, peg, new_egs, F)

        if F_p == [0,0]
            return th, new_egs
        end
    end

    return nothing, egs

end

function cp_not_move(c_p::myType.contact_pair, peg::myType.myPeg,
				egs::myType.mySocketEdges, config::myType.configSE2, F::Array, mode::Int)

	# rotate edge untill the contact points of CP not moving

    new_peg = Kinematics.get_peg_from_config(peg, config)
    n = length(egs.edges)
    if mode == 0
        px = new_peg.points[c_p.c_p.index].x
        py = new_peg.points[c_p.c_p.index].y
    elseif mode == 1
        if c_p.c_e.index > n/2
            px = egs.edges[c_p.c_e.index].e.x
            py = egs.edges[c_p.c_e.index].e.y
        else
            px = egs.edges[c_p.c_e.index].s.x
            py = egs.edges[c_p.c_e.index].s.y
        end
    end


    for th = 0:pi/90:pi

        new_egs = edge_rotate(egs, th, c_p.c_e.index, [px, py])
        #Draw.draw_open_sth(new_socket)
        F_p, Fn = Insertion.force_between_pair(c_p, peg, new_egs, F)

        if F_p == [0,0]
            return th, new_egs
        end
    end


    for th = 0:-pi/90:-pi

        new_egs = edge_rotate(egs, th, c_p.c_e.index, [px, py])
        #Draw.draw_open_sth(new_socket)
        F_p, Fn = Insertion.force_between_pair(c_p, peg, new_egs, F)

        if F_p == [0,0]
            return th, new_egs
        end
    end

    return nothing, egs

end


function cm_not_move(c_m::myType.contact_mode, peg::myType.myPeg,
				egs::myType.mySocketEdges, config::myType.configSE2,
				F::Array, mode::Int)

	# rotate edge untill the Contact Points of MoC not moving
    n = length(c_m.cps)

    th_list = Array{Array}(undef, 0)

    for i=1:n

        th, egs=cp_not_move(c_m.cps[i], peg, egs, config, F, mode)

        if th == nothing
            # println(i,"th cp has no th")
            continue
        end
        push!(th_list, [i, th])
    end

    return th_list, egs

end


function cm_not_move(c_m::myType.contact_mode, peg::myType.myPeg,
				socket::myType.mySocket, config::myType.configSE2, F::Array, mode::Int)

    n = length(c_m.cps)

    th_list = Array{Array}(undef, 0)

    for i=1:n

        th, socket=cp_not_move(c_m.cps[i], peg, socket, config, F, mode)

        if th == nothing
            # println(i,"th cp has no th")
            continue
        end
        push!(th_list, [i, th])
    end

    return th_list, socket

end

function edges_from_sink(socket::myType.mySocket, peg::myType.myPeg, F=[0,0,0],
	closed = false, convex = false, mode=1)

	# rotate edge until break the sink
	all_peg, c_p_s, model_list = Insertion.all_combinations(peg, socket, closed, convex)
	config_list =Insertion.cp_cf_list(peg, all_peg)
	cps, cms = Insertion.get_cp_cm(peg, socket, c_p_s, model_list)
	sink = Insertion.get_sink(peg, socket, F, closed, convex)

	n = length(sink)

	es = myType.mySocketEdges(socket)

	for i = 1:n
		m = sink[i]

		if F[1]==0 && F[2]==0
			F = Insertion.th_get_F(config_list[m].th, 2)
		end

		th_list, new_es = cm_not_move(cms[m], peg, es, config_list[m], F, mode)

		for j = 1:length(cms[m].cps)
			id = cms[m].cps[j].c_e.index
			es.edges[id] = new_es.edges[id]
		end
	end

	return es


end

function socket_from_edges(es::myType.mySocketEdges)
	n = length(es.edges)-1
	new_es = deepcopy(es)
	for i = 1:div(n,2)-1
		next = i+1
		x = 0
		y = 0
		if es.edges[i].e.y>es.edges[next].s.y
			l = get_line(es.edges[i].s, es.edges[i].e)
			x = es.edges[next].s.x
			if !all_close(l[2],0)
				y = (l[3]-l[1]*x)/l[2]
			else
				y = (es.edges[i].e.y+es.edges[i].s.y)/2
			end

		else
			l = get_line(es.edges[i].s, es.edges[i].e)
			y = es.edges[next].s.y
			if !all_close(l[1],0)
				x = (l[3]-l[2]*y)/l[1]
			else
				x = (es.edges[i].e.x+es.edges[i].s.x)/2
			end
		end
		v = [x-es.edges[next].s.x, y-es.edges[next].s.y]
		new_es.edges[i].e.x = x
		new_es.edges[i].e.y = y
		new_es.edges[next].s.x = x
		new_es.edges[next].s.y = y
		new_es.edges[next].e.x = es.edges[next].e.x + v[1]
		new_es.edges[next].e.y = es.edges[next].e.y + v[2]
	end

	for j = n:-1:n-div(n,2)+2
		last = j-1
		x = 0
		y = 0
		if es.edges[j].s.y>es.edges[last].e.y
			l = get_line(es.edges[j].s, es.edges[j].e)
			x = es.edges[last].e.x
			if !all_close(l[2],0)
				y = (l[3]-l[1]*x)/l[2]
			else
				y = (es.edges[j].e.y+es.edges[j].s.y)/2
			end

		else
			l = get_line(es.edges[j].s, es.edges[j].e)
			y = es.edges[last].e.y
			if !all_close(l[1],0)
				x = (l[3]-l[2]*y)/l[1]
			else
				x = (es.edges[j].e.x+es.edges[j].s.x)/2
			end
		end
		v = [x-es.edges[last].e.x, y-es.edges[last].e.y]
		new_es.edges[j].s.x = x
		new_es.edges[j].s.y = y
		new_es.edges[last].e.x = x
		new_es.edges[last].e.y = y
		new_es.edges[last].s.x = es.edges[last].s.x + v[1]
		new_es.edges[last].s.y = es.edges[last].s.y + v[2]
	end
	return new_es

end

function limit_th(peg::myType.myPeg, socket::myType.mySocket, F=[0,0,0],
	closed = false, convex = false)

	#limit rotate th for insertion

	es = myType.mySocketEdges(socket)
	n = length(es.edges)-1

	v1 = [es.edges[1].e.x - es.edges[1].s.x,
	es.edges[1].e.y - es.edges[1].s.y, 0]

	v2 = [es.edges[n].s.x - es.edges[n].e.x,
	es.edges[n].s.y - es.edges[n].e.y, 0]

	new_es = edges_from_sink(socket, peg, F, closed, convex)

	v3 = [new_es.edges[1].e.x - new_es.edges[1].s.x,
	new_es.edges[1].e.y - new_es.edges[1].s.y, 0]

	v4 = [new_es.edges[n].s.x - new_es.edges[n].e.x,
	new_es.edges[n].s.y - new_es.edges[n].e.y, 0]


	th1 = Insertion.compute_signed_angle(v1, v3)
	th2 = Insertion.compute_signed_angle(v4, v2)


	return min(th1, th2)


end

function limit_th_with_x_axis(peg::myType.myPeg, socket::myType.mySocket, F=[0,0,0],
	closed = false, convex = false)

	#limit rotate th for insertion


	new_es = edges_from_sink(socket, peg, F, closed, convex)
	n = length(new_es.edges)


	v1 = [new_es.edges[1].e.x - new_es.edges[1].s.x,
	new_es.edges[1].e.y - new_es.edges[1].s.y, 0]

	v2 = [new_es.edges[n].s.x - new_es.edges[n].e.x,
	new_es.edges[n].s.y - new_es.edges[n].e.y, 0]


	th1 = Insertion.compute_signed_angle(v1, [1, 0, 0])
	th2 = Insertion.compute_signed_angle([1, 0, 0], v2)


	return min(th1, th2)

end

function rotation_upper_limit(peg::myType.myPeg, socket::myType.mySocket, desired_sinks)
	fcone = friction_cone(socket)
	ths = CMT_graph.insertion_error(peg, socket, desired_sinks)
	sort!(ths)
	if length(ths) < 1
		return [-pi/2, pi/2]
	end
	th_min = ths[1]
	th_max = ths[length(ths)]
	d_fcone = abs((fcone[1][2] - fcone[1][1])/2)
	th_min = th_min - d_fcone
	th_max = th_max + d_fcone
	up_limit = [th_min + pi/2, th_max - pi/2]
	return up_limit
end

function entry_angle(socket)
	v1 = [socket.points[2].x - socket.points[1].x,
	socket.points[2].y - socket.points[1].y, 0]
	v2 = [socket.points[length(socket.points)-1].x - socket.points[length(socket.points)].x,
	socket.points[length(socket.points)-1].y - socket.points[length(socket.points)].y, 0]
	v3 = [0,-1,0]
	th1 = Insertion.compute_signed_angle(v3, v1)
	th2 = Insertion.compute_signed_angle(v3, v2)
	return [th1, th2]
end

function test_limit(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
	# cp_ids, cms, c_lists = ContactMode.all_combinations(peg,socket,c_ids)

	# println(cms)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
	up_limit = rotation_upper_limit(peg, socket, desired_sinks)
	ea = entry_angle(socket)
	# println(up_limit)
	# println(ea)
end

function eliminate_undesired_sink(peg, socket, c_ids, i = [], pic_count = 0, dt = 0.01)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)

	stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
	desired_sinks = CMT_graph.desired_sinks(c_ids, sinks)
	up_limit = rotation_upper_limit(peg, socket, desired_sinks)
	e_angle = entry_angle(socket)
	new_peg = deepcopy(peg)
	new_socket = deepcopy(socket)
	while stuck_sink != [] && abs(e_angle[1]-e_angle[2]) < abs(up_limit[1]-up_limit[2])
		t = Adjust.compute_t(new_peg, new_socket)
		new_socket = Rotation.adj_upsocket(new_socket, dt)
		new_peg = Adjust.get_peg_from_t(new_peg, new_socket, t)
		sinks = CMT_graph.get_sinks(new_peg,new_socket,c_ids)
		stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
		desired_sinks = CMT_graph.desired_sinks(c_ids, sinks)
		up_limit = rotation_upper_limit(new_peg, new_socket, desired_sinks)
		e_angle = entry_angle(new_socket)
		# println("i: ", i)
		Draw.draw_anime(new_peg, new_socket, i, pic_count)
		# PyPlot.clf()
		# PyPlot.axis("off")
		# Draw.draw_sth(new_peg, (1,0,0), 3.)
		# Draw.draw_open_sth(new_socket, (0.5,.4,.2), 3.)
		# CMT_graph.cmt_graph(new_peg,new_socket,c_ids)
		# println("udsinks: ", stuck_sink)
		pic_count += 1
	end

	if abs(e_angle[1]-e_angle[2]) >= abs(up_limit[1]-up_limit[2])
		println("can't eliminate undesired sinks")
		return nothing, nothing, nothing
	else
		return new_peg, new_socket, pic_count
	end
end

function test6(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
	new_peg, new_socket = eliminate_undesired_sink(peg, socket, c_ids, 0.1)
	sinks = CMT_graph.get_sinks(new_peg,new_socket,c_ids)
	println(sinks)
	Draw.draw_sth(new_peg)
	Draw.draw_open_sth(new_socket)
end

function test5(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
	# cp_ids, cms, c_lists = ContactMode.all_combinations(peg,socket,c_ids)

	# println(cms)
	t = Adjust.compute_t(peg, socket)
	socket = Rotation.adj_upsocket(socket, 0.6)
	peg = Adjust.get_peg_from_t(peg, socket, t)
	fcone = friction_cone(socket)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
	desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
	# println(stuck_sink)
	# println(desired_sinks)
	ths = CMT_graph.insertion_error(peg, socket, stuck_sink)
	# println("Thetas for sinks: ", ths)
	Fs = force_analysis_for_sink(socket, peg, stuck_sink[1], ths[1])
	# println("friction_cone: ", fcone)
	# println("Force of points: ", Fs)
	cm = myType.cm_id(stuck_sink[1])
	result, c_list= ContactMode.is_cm_valid(peg, socket, cm)
	new_peg = Kinematics.get_peg_from_config(peg, c_list[1])
	Draw.draw_sth(new_peg)
	Draw.draw_open_sth(socket)
end

function test4(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	closed = false
	convex = Insertion.is_convex(socket)
	F=[0,0,0]
	es = edges_from_sink(socket, peg, F, closed, convex)
	th = limit_th(peg, socket, F, closed, convex)
	# println(th)
	#new_es = socket_from_edges(es)
	# Draw.draw_open_sth(socket)
	# Draw.draw_edges(new_es, (0,0,1), 2.0)
	# PyPlot.figure()
	Draw.draw_open_sth(socket)
	Draw.draw_edges(es, (0,0,1), 2.0)
end



function test1(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    Draw.draw_open_sth(socket)
    pc = [socket.points[2].x, socket.points[2].y]

    egs = myType.mySocketEdges(socket)
    new_egs = edge_rotate(egs, pi/6, 1, pc)

    Draw.draw_edges(new_egs, (0,0,1), 2.0)

end

function test2(i)
        peg, socket = Construct_Peg_Socket.read_from_file(i)

    	all_peg, c_p_s, model_list = Insertion.all_combinations(peg, socket)

    	#Draw.draw_all_peg(all_peg, peg, socket)

    	config_list =Insertion.cp_cf_list(peg, all_peg)

    	# println(c_p_s)
    	# println(model_list)
        egs = myType.mySocketEdges(socket)

    	cps, cms = Insertion.get_cp_cm(peg, socket, c_p_s, model_list)

        F = Insertion.th_get_F(config_list[3].th, 2)
        mode = 1
        th, new_egs=cp_not_move(cms[3].cps[1], peg, egs, config_list[3], F, mode)

        if th == nothing
            println("no th")
        end
        Draw.draw_open_sth(socket)
        Draw.draw_edges(new_egs, (0,0,1), 2.0)

end

function test3(i)
        peg, socket = Construct_Peg_Socket.read_from_file(i)

    	all_peg, c_p_s, model_list = Insertion.all_combinations(peg, socket)

    	#Draw.draw_all_peg(all_peg, peg, socket)

    	config_list =Insertion.cp_cf_list(peg, all_peg)

    	# println(c_p_s)
    	# println(model_list)

        egs = myType.mySocketEdges(socket)

    	cps, cms = Insertion.get_cp_cm(peg, socket, c_p_s, model_list)

        F = Insertion.th_get_F(config_list[29].th, 2)
        mode = 1
        th_list, new_egs=cm_not_move(cms[29], peg, egs, config_list[29], F, mode)

        println(th_list)
        PyPlot.figure()
        Draw.draw_open_sth(socket)
        Draw.draw_edges(new_egs, (0,0,1), 2.0)

end

end
