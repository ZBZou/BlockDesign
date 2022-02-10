module CMT_graph

using ..myType
using Compose
using Gadfly
using Colors
using GR
using Combinatorics

import LightGraphs
import GraphPlot
import GraphRecipes
import LinearAlgebra
import JuMP
import Cairo
import Fontconfig

import ..CONST
import ..Construct_Peg_Socket
import ..Constraints
import ..Kinematics
import ..ConcaveTest
import ..Draw
import ..ContactMode
import ..Insertion
import ..Adjust
import ..Rotation

import PyPlot
import Plots
import GraphPlot

function all_close(a, b, ep=0.0001)
	if abs(a - b) < ep
		return true
	end
	return false;
end


function CM_tansition_with_one_CP(cm1::myType.contact_mode, cm2::myType.contact_mode)
    if !ContactMode.is_pp_contact(cm2)
        return true
    else
        if !inverse_backchining_sequence(cm1, cm2)
            return false
        else
            return true
        end
    end
end

function inverse_backchining_sequence(cm1::myType.contact_mode, cm2::myType.contact_mode)
    v1 = [0,-1]
    e_ps = [cm1.cps[1].c_e.s, cm1.cps[1].c_e.e]

    e_p_ids = [cm2.cps[1].c_e.s.index, cm2.cps[1].c_e.e.index,
    cm2.cps[2].c_e.s.index, cm2.cps[2].c_e.e.index]

    id = 0
    p1 = 0
    p2 = 0
    for i = 1:length(e_p_ids)
        if i == 1
            id = e_p_ids[i]
            continue
        end

        if e_p_ids[i] != id
            id = e_p_ids[i]
            continue
        else
            break
        end
    end

    for p in e_ps
        if p.index == id
            p2 = p
        else
            p1 = p
        end
    end

    v2 = [p2.x - p1.x, p2.y-p1.y]

    if LinearAlgebra.dot(v1, v2) >=0
        return true
    else
        return false
    end

end

function is_neighbor_modes(cm1::myType.contact_mode, cm2::myType.contact_mode)

	# test if two modes are neighbors

	 l1 = length(cm1.cps)
	 l2 = length(cm2.cps)
	#
	# if abs(l1 - l2) != 1
	# 	return false;
	# end

	cps1 = Array{Int}(undef, 0)
	cps2 = Array{Int}(undef, 0)

	sum = Array{Int}(undef, 0)

	for i = 1:l1
		push!(cps1, cm1.cps[i].id)
		push!(sum, cm1.cps[i].id)
	end
	for i = 1:l2
		push!(cps2, cm2.cps[i].id)
		push!(sum, cm2.cps[i].id)
	end

	r1 = false
	r2 = false

	for i = 1:l1
		if !(cps1[i] in cps2)
			r1 = true
		end
	end

	for i = 1:l2
		if !(cps2[i] in cps1)
			r2 = true
		end
	end

	if r1 && r2 || !r1 && !r2
		return false
	end


	inter = intersect(cps1, cps2)
	# println("here1")
	# println(cps1)
	# println(cps2)
	 # println(inter)
	 # println(sum)

	l3 = length(sum)

	for i = l3:-1:1

		if sum[i] in inter
			deleteat!(sum, i)
		end
	end
	# println("here3")
	# println(sum)

	if isempty(sum)
		return false
	end


	set_p = Array{Int}(undef, 0)
	for i = 1: length(sum)
		if r1
			for j = 1:length(cm1.cps)
				if cm1.cps[j].id == sum[i]
					push!(set_p, cm1.cps[j].c_p.index)
				end
			end
		elseif r2
			for j = 1:length(cm2.cps)
				if cm2.cps[j].id == sum[i]
					push!(set_p, cm2.cps[j].c_p.index)
				end
			end
		else
			prinln("error in neighbor")
		end
	end

	set_p = unique!(set_p)



	if length(set_p) != 1
		return false;
	end
	# else if the number of difference is larger than 1
	# the two modes cannot be neighbors

	# sort_contact_mode(cm1)
	# sort_contact_mode(cm2)

	# after sort
	# compare

	return true

end

function find_transition_configs(c_lists1::Array, c_lists2::Array, th_epsilon = 0.001)
    c1_min = c_lists1[1].th
    c1_max = c_lists1[length(c_lists1)].th

    c2_min = c_lists2[1].th
    c2_max = c_lists2[length(c_lists2)].th

    if c1_max < c2_min
        if c1_max + th_epsilon> c2_min
            return c_lists1[length(c_lists1)]
        else
            return nothing
        end
    elseif c1_min > c2_max
        if c1_min - th_epsilon < c2_max
            return c_lists1[1]
        else
            return nothing
        end
    end

    th_min = max(c1_min, c2_min)
    th_max = min(c1_max, c2_max)

    for config in c_lists1
        if all_close(config.th, (th_min+th_max)/2, th_epsilon)
            return config
        end
    end

end
#

function get_cp_cm(peg::myType.myPeg, socket::myType.mySocket,
	c_p_s::Array, model_list::Array)

	cps = Array{myType.contact_pair}(undef, 0)
	cms = Array{myType.contact_mode}(undef, 0)
	ces = myType.mySocketEdges(socket)





	for i = 1:length(c_p_s)
		cp_s = Array{myType.contact_pair}(undef, 0)
		push!(cps, myType.contact_pair(i, peg.points[c_p_s[i][1]], ces.edges[c_p_s[i][2]]))
	end



	for i = 1:length(model_list)
		cp_s = Array{myType.contact_pair}(undef, 0)
		for k = 1:length(model_list[i])
			push!(cp_s, cps[model_list[i][k]])
		end
		push!(cms, myType.contact_mode(cp_s))
	end

	return cps, cms

end
#
# function model_graph(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{myType.contact_id}, F::Array,
#     closed = false, convex = false)
#
# 	#find insertion CMT graph
#
#
# 	f_with_p = false
#
# 	if F[1]==0&&F[2]==0
# 		f_with_p = true
# 	end
#
#
#
#
# 	c_p_s, model_list, c_lists = ContactMode.all_combinations(peg, socket, c_ids)
#
#
# 	println(c_p_s)
# 	println(model_list)
#
#
# 	cps, cms = get_cp_cm(peg, socket, c_p_s, model_list)
# 	n = length(cms)
# 	g = LightGraphs.DiGraph(n+1)
#
#     # println(cms)
#     println(n)
#
# 	# if valid_direction(cms[6], cms[11], peg, socket, config_list[6], F)
# 	# 	println("here")
# 	# end
#
# 	cmp_id = Array{Array}(undef, 0)
#
# 	for i = 1:n
# 		for j = 1:n
# 			#println([i,j])
#
# 			result, id = Insertion.is_p_p(cms[i])
# 			if result
# 				push!(cmp_id, [i,id])
# 				continue
# 			end
#
#             config = find_transition_configs(c_lists[i], c_lists[j])
#
# 			if i!=j && Insertion.is_neighbor_modes(cms[i], cms[j]) && config!=nothing
#
# 				if f_with_p
# 				 F = Insertion.th_get_F(config.th, 2)
# 			 	end
# 				# println(F)
# 				# println(valid_direction(cms[i], cms[j], peg, socket, config_list[i], F))
# 				if Insertion.valid_direction(cms[i], cms[j], peg, socket, config, F)
#
# 					LightGraphs.add_edge!(g, i, j)
# 				end
# 			end
# 		end
# 	end
#
#
# 	for i = 1:length(cmp_id)
# 		for j = 1: n
#             config = find_transition_configs(c_lists[cmp_id[i][1]], c_lists[j])
#
# 			if i!=j && Insertion.is_neighbor_modes(cms[cmp_id[i][1]], cms[j])&& config!=nothing
#
# 				if f_with_p
# 					F = Insertion.th_get_F(config.th, 2)
# 				end
#
# 				if Insertion.is_p_p_pair(cms[cmp_id[i][1]], cms[j])
#
# 					if Insertion.convex_pp(socket, cms[cmp_id[i][1]], cms[j])
# 						result1,result2 = Insertion.p_p_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config, F, cmp_id[i][2])
# 					else
# 						result1,result2 = Insertion.p_p_concave(cms[cmp_id[i][1]], cms[j], peg, socket, config, F, cmp_id[i][2])
# 					end
#
# 					if result1
#
# 						LightGraphs.add_edge!(g, cmp_id[i][1], j)
# 					end
#
# 					if result2
#
# 						LightGraphs.add_edge!(g, cmp_id[i][1], n+1)
#
# 					end
# 				elseif Insertion.valid_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config, F)
#
# 					LightGraphs.add_edge!(g, cmp_id[i][1], j)
# 				end
#
# 			end
# 		end
#
# 	end
#
#
# 	if f_with_p
# 		for i = 1:length(c_p_s)
# 			LightGraphs.add_edge!(g, n+1, i)
# 		end
#
# 	else
# 		id_set = Insertion.cps_trans_id(socket, F, c_p_s)
#
# 		for i = 1:length(id_set)
# 			LightGraphs.add_edge!(g, n+1, i)
# 		end
# 	end
#
# 	# for i = 1:length(id_set)
# 	# 	LightGraphs.add_edge!(g, n+1, i)
# 	# end
#
#
# 	for i =1:length(cps)
# 		if f_with_p
#             config = c_lists[i][floor(Int,(1+length(c_lists))/2)]
# 			F = Insertion.th_get_F(config.th, 2)
# 		end
# 		F_p, Fn = Insertion.force_between_pair(cps[i], peg, socket, config, F)
#
# 		if Fn<0
# 			LightGraphs.add_edge!(g, i, n+1)
# 		end
# 	end
#
#
#
# 	LightGraphs.collect(LightGraphs.edges(g))
# 	A = LightGraphs.LinAlg.adjacency_matrix(g)
#
#
# 	s_A = size(A, 1)
# 	sink = Array{Int}(undef, 0)
# 	for i = 1:s_A
# 		num = 0
# 		if !(i in A.rowval)
# 			push!(sink, i)
# 		end
# 	end
# 	println("sink:", sink)
#
#
#
#
# 	# if count == 0
# 	# 	error("no edge")
# 	#
# 	# else
#
# 	# GraphRecipes.graphplot(g, names=1:n+1,  linecolor = :darkgrey, color = :white,
# 	# 							markersize=1.5, arrow=Plots.arrow(:closed, :head, 1, 1), shorten=0.1)
#
# 	nodesize = Vector{Float64}(undef, 0)
# 	nodelabelsize = Vector{Float64}(undef, 0)
# 	membership = Array{Int}(undef, 0)
# 	nodecolor = [colorant"lightseagreen", colorant"orange"]
# 	edgestrokec = colorant"black"
#
# 	for i = 1:n+1
# 		push!(nodesize, 1.)
# 		push!(nodelabelsize, 1.)
# 		push!(membership, 1)
# 	end
#
# 	for j in sink
# 		membership[j] = 2
# 	end
# 	nodefillc = nodecolor[membership]
#
# 	GraphPlot.draw(PNG("insertion_graph.png", 16cm, 16cm), GraphPlot.gplot(g, nodelabel=1:n+1,
# 	arrowlengthfrac=0.05,
# 	nodesize=nodesize, nodelabelsize=nodelabelsize, nodefillc=nodefillc,
# 	edgestrokec=edgestrokec))
# 	return g
#
# end


function force_analysis(peg, socket, cm, config, F)

    app = []
    lea = []
    for k = 1:length(cm.cps)
        push!(app, false)
        push!(lea, false)
    end

    for i =1:length(cm.cps)
        target_p = cm.cps[i].c_p
        F_all = deepcopy(F)
        stuck = false
        for j = 1:length(cm.cps)
			if target_p.index == cm.cps[j].c_p.index
				continue
			end

			F_j = deepcopy(F)

			F_pj, Fnj= Insertion.force_between_pair(cm.cps[j], peg, socket, F_j)

            if F_pj[1] == 0 && F_pj[2] == 0
                stuck = true
                break
            end

			if Fnj>0

				F_all[1] += F_pj[1]
                F_all[2] += F_pj[2]
			end
		end

        if stuck == true
            continue
        end

        F_p, Fn = Insertion.force_between_pair(cm.cps[i], peg, socket, F_all)

		if Fn < 0
            lea[i] = true
        else
            app[i] = true
        end

    end

    return app, lea

end



function corner_out(peg, socket, cm1, cm2, config, F, p_id)

    F_all = deepcopy(F)
    cp_id = 0
	r = false
	es = myType.mySocketEdges(socket)
	for i = 1:length(cm2.cps)
		if cm2.cps[i].c_p.index == p_id
			cp_id = i
			continue
		else
			r = true
		end
	end

	if !r
		return false
	end


    for j = 1:length(cm2.cps)

		if cm2.cps[j].c_p.index == p_id

			continue
		end

		F_j = deepcopy(F)

		F_pj, Fnj= Insertion.force_between_pair(cm2.cps[j], peg, socket, F_j)

        if F_pj[1] == 0 && F_pj[2] == 0
			println("here1")
            return false
        end

		if Fnj>0

			F_all[1] += F_pj[1]
            F_all[2] += F_pj[2]
		end
	end

	e_id = []
	for i = 1:length(cm1.cps)
		if cm1.cps[i].c_p.index == p_id
			push!(e_id, cm1.cps[i].c_e.index)
		end
	end

	sort!(e_id)
	v1 = [es.edges[e_id[1]].s.x - es.edges[e_id[1]].e.x,
	es.edges[e_id[1]].s.y - es.edges[e_id[1]].e.y, 0]

	v2 = [es.edges[e_id[2]].e.x - es.edges[e_id[2]].s.x,
	es.edges[e_id[2]].e.y - es.edges[e_id[2]].s.y, 0]

	v = [0,0,0]

	for i = 1:length(cm2.cps)
		if cm2.cps[i].c_p.index == p_id
			if e_id[1] == cm2.cps[i].c_e.index
				v = v1
			elseif e_id[2] == cm2.cps[i].c_e.index
				v = v2
			end
		end
	end

	if LinearAlgebra.dot(v, F_all) < 0
		# println("here2")
		# println([v, F_all])
		return false
	end
	if cp_id == 0
		return false
	end
    F_p, Fn = Insertion.force_between_pair(cm2.cps[cp_id], peg, socket, F_all)

	if Fn < 0 || (F_p[1] == 0 && F_p[2] == 0)
		# println("here3")
        return false
    end

    return true
end


function corner_in(peg, socket, cm1, cm2, config, F, p_id)
    F_all = deepcopy(F)
    cp_id = 0
	r = false
	es = mySocketEdges(socket)
	for i = 1:length(cm2.cps)
		if cm2.cps[i].c_p.index == p_id
			cp_id = i
			continue
		else
			r = true
		end
	end

	if !r
		return true
	end


    for j = 1:length(cm2.cps)

		if cm2.cps[j].c_p.index == p_id

			continue
		end

		F_j = deepcopy(F)

		F_pj, Fnj= Insertion.force_between_pair(cm2.cps[j], peg, socket, F_j)

        if F_pj[1] == 0 && F_pj[2] == 0

            return false
        end

		if Fnj>0

			F_all[1] += F_pj[1]
            F_all[2] += F_pj[2]
		end
	end

	e_id = []
	for i = 1:length(cm1.cps)
		if cm1.cps[i].c_p.index == p_id
			push!(e_id, cm1.cps[i].c_e.index)
		end
	end

	sort!(e_id)
	v1 = [es.edges[e_id[1]].s.x - es.edges[e_id[1]].e.x,
	es.edges[e_id[1]].s.y - es.edges[e_id[1]].e.y, 0]

	v2 = [es.edges[e_id[2]].e.x - es.edges[e_id[2]].s.x,
	es.edges[e_id[2]].e.y - es.edges[e_id[2]].s.y, 0]

	v = [0,0,0]

	for i = 1:length(cm2.cps)
		if cm2.cps[i].c_p.index == p_id
			if e_id[1] == cm2.cps[i].c_e.index
				v = v1
			elseif e_id[2] == cm2.cps[i].c_e.index
				v = v2
			end
		end
	end

	if LinearAlgebra.dot(v, F_all) > 0
		return false
	end
	if length(cm2.cps) < cp_id
		return false
	end
	if cp_id == 0
		return false
	end
    F_p, Fn = Insertion.force_between_pair(cm2.cps[cp_id], peg, socket, F_all)

	if Fn < 0 || (F_p[1] == 0 && F_p[2] == 0)
        return false
    end

    return true
end

function get_cp_cm(peg::myType.myPeg, socket::myType.mySocket,
	c_p_s::Array, model_list::Array)

	cps = Array{myType.contact_pair}(undef, 0)
	cms = Array{myType.contact_mode}(undef, 0)
	ces = myType.mySocketEdges(socket)





	for i = 1:length(c_p_s)
		cp_s = Array{myType.contact_pair}(undef, 0)
		push!(cps, myType.contact_pair(i, peg.points[c_p_s[i].p_id], ces.edges[c_p_s[i].e_id]))

	end



	for i = 1:length(model_list)


		cp_s = Array{myType.contact_pair}(undef, 0)
		for k = 1:length(model_list[i])
			push!(cp_s, cps[model_list[i][k]])
		end
		push!(cms, myType.contact_mode(cp_s))



	end

	return cps, cms


end


function is_p_p(cm1)
	# test if point-point situation

	n = length(cm1.cps)
	id = 0
	for i = 1:n
		for j = 1:n
			if i!=j && cm1.cps[i].c_p.index == cm1.cps[j].c_p.index
				id = cm1.cps[i].c_p.index
				return true,id
			end
		end
	end
	return false,id
end

function p_p_id(cm::myType.contact_mode)
	temp = []
	p_id = 0
	for i = 1:length(cm.cps)
		if !(cm.cps[i].c_p.index in temp)
			push!(temp, cm.cps[i].c_p.index)
		else
			p_id = cm.cps[i].c_p.index
		end
	end
	return p_id
end


function convex_pp(socket::myType.mySocket, cm::myType.contact_mode)
	es = myType.mySocketEdges(socket)
	temp = []
	p_id = 0
	for i = 1:length(cm.cps)
		if !(cm.cps[i].c_p.index in temp)
			push!(temp, cm.cps[i].c_p.index)
		else
			p_id = cm.cps[i].c_p.index
		end
	end

	e_id = []
	for i = 1:length(cm.cps)
		if cm.cps[i].c_p.index == p_id
			push!(e_id, cm.cps[i].c_e.index)
		end
	end

	sort!(e_id)
	v1 = [es.edges[e_id[1]].s.x - es.edges[e_id[1]].e.x,
	es.edges[e_id[1]].s.y - es.edges[e_id[1]].e.y, 0]

	v2 = [es.edges[e_id[2]].s.x - es.edges[e_id[2]].e.x,
	es.edges[e_id[2]].s.y - es.edges[e_id[2]].e.y, 0]

	th = compute_signed_angle(v1, v2)

	if th > 0
		return true
	else
		return false
	end
end


function compute_signed_angle(v1, v2)

	vn = [0,0,1]


	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return th


end

function cm_to_contact_ids(cm::myType.contact_mode)
	c_ids = []
	n = length(cm.cps)
	for i = 1:n
		c_id = myType.contact_id(cm.cps[i].c_p.index, cm.cps[i].c_e.index)
		push!(c_ids, c_id)
	end
	return c_ids
end

function cps_to_contact_ids(lists, cps::Array)
	c_ids = []
	for id in lists
		push!(c_ids, cps[id])
	end
	return c_ids
end

function outward_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(-pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

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

function force_analysis_for_sink(socket, peg, sink::Array, theta)
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

function two_p_one_e(cm::Array)
	n = length(cm)
	ne = []
	for i = 1:n
		push!(ne, cm[i].e_id)
	end
	unique!(ne)
	if length(ne) < n
		return true
	else
		return false
	end
end

function two_p_one_e(cm::contact_mode)
	n = length(cm.cps)
	ne = []
	for i = 1:n
		push!(ne, cm.cps[i].c_e.index)
	end
	unique!(ne)
	if length(ne) < n
		return true
	else
		return false
	end
end

function rotate_sink(cm)
	body
end

function is_in_friction_cone(peg, socket, sink::Array)
	fcone = friction_cone(socket)
	cm = myType.cm_id(sink)
	th = insertion_error(peg, socket, cm)
	Fs = force_analysis_for_sink(socket, peg, sink, th)
	n = length(sink)
	p_id = []
	for i = 1:n
		f_i = Fs[i]
		fcone_i = fcone[sink[i].e_id]
		if f_i >= fcone_i[1] && f_i <= fcone_i[2]
			push!(p_id, sink[i].p_id)
		end
	end
	return p_id
end

function minimul_cp_distance(c_ids, sink)

end


function rotate_direction(args)
	p_c = Insertion.compute_centroid(peg)
end

function test5(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
	# cp_ids, cms, c_lists = ContactMode.all_combinations(peg,socket,c_ids)

	# println(cms)
	fcone = friction_cone(socket)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
	desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
	# println(stuck_sink[1])
	p_id = is_in_friction_cone(peg, socket, stuck_sink[1])

	println(stuck_sink)
	println(desired_sinks)
	ths = CMT_graph.insertion_error(peg, socket, stuck_sink)
	println("Thetas for sinks: ", ths)
	Fs = force_analysis_for_sink(socket, peg, stuck_sink[1], ths[1])
	println("friction_cone: ", fcone)
	println("Force of points: ", Fs)
	println("p_id: ", p_id)
	cm = myType.cm_id(stuck_sink[1])
	result, c_list= ContactMode.is_cm_valid(peg, socket, cm)
	new_peg = Kinematics.get_peg_from_config(peg, c_list[1])
	Draw.draw_sth(new_peg)
	Draw.draw_open_sth(socket)
end


function cmt_graph(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{myType.contact_id}, F=[0,0,0],
    closed = false, convex = false)

    c_p_s, cm_ids, c_lists = ContactMode.find_valid_three_pairs_cms(peg, socket, c_ids)

    all_cms = ContactMode.contact_modes_from_three_pairss(cm_ids)
    cps, cms = get_cp_cm(peg, socket, c_p_s, all_cms)

    l_all = length(all_cms)
	g = LightGraphs.DiGraph(l_all+1)
    n = length(cm_ids)

    f_with_p = false

    if F[1] == 0 && F[2] == 0
        f_with_p = true
    end

    for i = 1:n

        config = c_lists[i][1]
        # cms_i = []
        if f_with_p
            F = Insertion.th_get_F(config.th, 2)
        end

        temp1 = collect(combinations(cm_ids[i],1))
        temp2 = collect(combinations(cm_ids[i],2))

        id_i = 0

        for m = 1:l_all
            if cm_ids[i] == all_cms[m]
                id_i = m
            end
        end

        result,pp_id = is_p_p(cms[id_i])
        if !result
			r1 = false
			cm_i = cm_to_contact_ids(cms[id_i])
			r1 = two_p_one_e(cm_i)

			p_id = is_in_friction_cone(peg, socket, cm_i)
            app,lea = force_analysis(peg, socket, cms[id_i], config, F)
            # println(lea)
            for j = 1:length(temp2)
				cm_j = cps_to_contact_ids(temp2[j], c_p_s)
                id_j = 0
                for m = 1:l_all
                    if temp2[j] == all_cms[m]
                        id_j = m
                    end
                end

                for k = 1:length(cm_ids[i])
                    if !(k in temp2)
                        if app[k] == true
                            LightGraphs.add_edge!(g, id_j, id_i)
                        end
                        if lea[k] == true || (r1 && length(p_id) == 0 && is_goal_cm(c_ids, cm_j))
                            LightGraphs.add_edge!(g, id_i, id_j)
                        end
                    end
                end

                app_j,lea_j = force_analysis(peg, socket, cms[id_j], config, F)

                for k = 1:length(temp2[j])
                    id_k = 0
                    for m = 1:l_all
                        if [temp2[j][k]] == all_cms[m]
                            id_k = m
                        end
                    end
                    LightGraphs.add_edge!(g, l_all+1, id_k)
                    LightGraphs.add_edge!(g, id_k, id_j)
                    if lea_j[k] == true
                        LightGraphs.add_edge!(g, id_j, id_k)
                    end
					if r1 && (length(p_id) == 1) && (p_id[1] == temp2[j][k].p_id)
						LightGraphs.add_edge!(g, id_i, id_k)
					end
                end
            end
        else
            if !convex_pp(socket, cms[id_i])
					# app,lea = force_analysis(peg, socket, cms[id_i], config, F)
	            for j = 1:length(temp2)
	                id_j = 0
	                for m = 1:l_all
	                    if temp2[j] == all_cms[m]
	                        id_j = m
	                    end
	                end

	                LightGraphs.add_edge!(g, id_j, id_i)
	                result2, id_nothing = Insertion.is_p_p(cms[id_j])

	                if result2
	                    LightGraphs.add_edge!(g, id_j, l_all+1)
	                    LightGraphs.add_edge!(g, l_all+1, id_j)
	                    for k = 1:length(temp2[j])
	                        id_k = 0
	                        for m = 1:l_all
	                            if [temp2[j][k]] == all_cms[m]
	                                id_k = m
	                            end
	                        end
	                        LightGraphs.add_edge!(g, l_all+1, id_k)
	                        LightGraphs.add_edge!(g, id_k, id_j)
	                        LightGraphs.add_edge!(g, id_j, id_k)
	                    end
	                else

	                    LightGraphs.add_edge!(g, id_i, id_j)

	                    for k = 1:length(temp2[j])
	                        id_k = 0
	                        for m = 1:l_all
	                            if [temp2[j][k]] == all_cms[m]
	                                id_k = m
	                            end
	                        end
	                        LightGraphs.add_edge!(g, id_k, id_j)
	                        LightGraphs.add_edge!(g, l_all+1, id_k)

	                        if cms[id_k].cps[1].c_p.index != pp_id
	                            LightGraphs.add_edge!(g, id_i, id_k)
	                            LightGraphs.add_edge!(g, id_j, id_k)
	                        end
	                    end
	                end
				end

			else
				pp_id = p_p_id(cms[id_i])

				for j = 1:length(temp2)

	                id_j = 0
	                for m = 1:l_all
	                    if temp2[j] == all_cms[m]
	                        id_j = m
	                    end
	                end


					if corner_in(peg, socket, cms[id_i], cms[id_j], config, F, pp_id)
						LightGraphs.add_edge!(g, id_j, id_i)
					end
					# println(cm_ids[i])
					if corner_out(peg, socket, cms[id_i], cms[id_j], config, F, pp_id)
						LightGraphs.add_edge!(g, id_i, id_j)
					end



					app_j,lea_j = force_analysis(peg, socket, cms[id_j], config, F)
	                for k = 1:length(temp2[j])
	                    id_k = 0
	                    for m = 1:l_all
	                        if [temp2[j][k]] == all_cms[m]
	                            id_k = m
	                        end
	                    end
	                    LightGraphs.add_edge!(g, l_all+1, id_k)
	                    LightGraphs.add_edge!(g, id_k, id_j)
	                    if lea_j[k] == true
	                        LightGraphs.add_edge!(g, id_j, id_k)
	                    end
	                end
				end

            end

        end



    end



	LightGraphs.collect(LightGraphs.edges(g))
	A = LightGraphs.LinAlg.adjacency_matrix(g)


	s_A = size(A, 1)
	sink = Array{Int}(undef, 0)
    sink_model = []
    pegs = []
    # for i = 1:n
    #     new_peg = Kinematics.get_peg_from_config(peg, c_lists[i][1])
    #     push!(pegs, new_peg)
    # end

	for i = 1:s_A
		num = 0
		if !(i in A.rowval)
			push!(sink, i)
            push!(sink_model, all_cms[i])
            for j = 1:n
                if all_cms[i] == cm_ids[j]
                    new_peg = Kinematics.get_peg_from_config(peg, c_lists[j][1])
                    push!(pegs, new_peg)
                end
            end
		end
	end

	# println("sink:", sink)
    # println("sink_model:", sink_model)

    Draw.draw_all_peg(pegs, socket)


	# if count == 0
	# 	error("no edge")
	#
	# else

	# GraphRecipes.graphplot(g, names=1:n+1,  linecolor = :darkgrey, color = :white,
	# 							markersize=1.5, arrow=Plots.arrow(:closed, :head, 1, 1), shorten=0.1)

	nodesize = Vector{Float64}(undef, 0)
	nodelabelsize = Vector{Float64}(undef, 0)
	membership = Array{Int}(undef, 0)
	nodecolor = [colorant"lightseagreen", colorant"orange"]
	edgestrokec = colorant"black"

	for i = 1:l_all+1
		push!(nodesize, 1.)
		push!(nodelabelsize, 1.)
		push!(membership, 1)
	end

	for j in sink
		membership[j] = 2
	end

	nodefillc = nodecolor[membership]

	GraphPlot.draw(PNG("./insertion_graph.png", 16cm, 16cm), GraphPlot.gplot(g, nodelabel=1:l_all+1,
	arrowlengthfrac=0.05,
	nodesize=nodesize, nodelabelsize=nodelabelsize, nodefillc=nodefillc,
	edgestrokec=edgestrokec))
	return g
end

function get_sinks(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{myType.contact_id}, F=[0,0,0],
    closed = false, convex = false)
    c_p_s, cm_ids, c_lists = ContactMode.find_valid_three_pairs_cms(peg, socket, c_ids)

    all_cms = ContactMode.contact_modes_from_three_pairss(cm_ids)
    cps, cms = get_cp_cm(peg, socket, c_p_s, all_cms)

    l_all = length(all_cms)
	g = LightGraphs.DiGraph(l_all+1)
    n = length(cm_ids)

    f_with_p = false

    if F[1] == 0 && F[2] == 0
        f_with_p = true
    end
	

    for i = 1:n

        config = c_lists[i][1]
        # cms_i = []
        if f_with_p
            F = Insertion.th_get_F(config.th, 2)
        end

        temp1 = collect(combinations(cm_ids[i],1))
        temp2 = collect(combinations(cm_ids[i],2))

        id_i = 0

        for m = 1:l_all
            if cm_ids[i] == all_cms[m]
                id_i = m
            end
        end

        result,pp_id = is_p_p(cms[id_i])
        if !result
			cm_i = cm_to_contact_ids(cms[id_i])
			r1 = two_p_one_e(cm_i)
			p_id = is_in_friction_cone(peg, socket, cm_i)
            app,lea = force_analysis(peg, socket, cms[id_i], config, F)
            # println(lea)
			# println(r1)
            for j = 1:length(temp2)
				cm_j = cps_to_contact_ids(temp2[j], c_p_s)
                id_j = 0
                for m = 1:l_all
                    if temp2[j] == all_cms[m]
                        id_j = m
                    end
                end

                for k = 1:length(cm_ids[i])
                    if !(k in temp2)
                        if app[k] == true
                            LightGraphs.add_edge!(g, id_j, id_i)
                        end
                        if lea[k] == true || (r1 && length(p_id) == 0 && is_goal_cm(c_ids, cm_j))
                            LightGraphs.add_edge!(g, id_i, id_j)
                        end
                    end
                end

                app_j,lea_j = force_analysis(peg, socket, cms[id_j], config, F)

                for k = 1:length(temp2[j])
                    id_k = 0
                    for m = 1:l_all
                        if [temp2[j][k]] == all_cms[m]
                            id_k = m
                        end
                    end
                    LightGraphs.add_edge!(g, l_all+1, id_k)
                    LightGraphs.add_edge!(g, id_k, id_j)
                    if lea_j[k] == true
                        LightGraphs.add_edge!(g, id_j, id_k)
                    end
					if r1 && length(p_id) == 1 && p_id[1] == temp2[j][k].p_id
						LightGraphs.add_edge!(g, id_i, id_k)
					end
                end
            end
        else
            if !convex_pp(socket, cms[id_i])
					# app,lea = force_analysis(peg, socket, cms[id_i], config, F)
	            for j = 1:length(temp2)
	                id_j = 0
	                for m = 1:l_all
	                    if temp2[j] == all_cms[m]
	                        id_j = m
	                    end
	                end

	                LightGraphs.add_edge!(g, id_j, id_i)
	                result2, id_nothing = Insertion.is_p_p(cms[id_j])

	                if result2
	                    LightGraphs.add_edge!(g, id_j, l_all+1)
	                    LightGraphs.add_edge!(g, l_all+1, id_j)
	                    for k = 1:length(temp2[j])
	                        id_k = 0
	                        for m = 1:l_all
	                            if [temp2[j][k]] == all_cms[m]
	                                id_k = m
	                            end
	                        end
	                        LightGraphs.add_edge!(g, l_all+1, id_k)
	                        LightGraphs.add_edge!(g, id_k, id_j)
	                        LightGraphs.add_edge!(g, id_j, id_k)
	                    end
	                else

	                    LightGraphs.add_edge!(g, id_i, id_j)

	                    for k = 1:length(temp2[j])
	                        id_k = 0
	                        for m = 1:l_all
	                            if [temp2[j][k]] == all_cms[m]
	                                id_k = m
	                            end
	                        end
	                        LightGraphs.add_edge!(g, id_k, id_j)
	                        LightGraphs.add_edge!(g, l_all+1, id_k)

	                        if cms[id_k].cps[1].c_p.index != pp_id
	                            LightGraphs.add_edge!(g, id_i, id_k)
	                            LightGraphs.add_edge!(g, id_j, id_k)
	                        end
	                    end
	                end
				end

			else
				pp_id = p_p_id(cms[id_i])

				for j = 1:length(temp2)

	                id_j = 0
	                for m = 1:l_all
	                    if temp2[j] == all_cms[m]
	                        id_j = m
	                    end
	                end


					if corner_in(peg, socket, cms[id_i], cms[id_j], config, F, pp_id)
						LightGraphs.add_edge!(g, id_j, id_i)
					end
					# println(cm_ids[i])
					if corner_out(peg, socket, cms[id_i], cms[id_j], config, F, pp_id)
						LightGraphs.add_edge!(g, id_i, id_j)
					end



					app_j,lea_j = force_analysis(peg, socket, cms[id_j], config, F)
	                for k = 1:length(temp2[j])
	                    id_k = 0
	                    for m = 1:l_all
	                        if [temp2[j][k]] == all_cms[m]
	                            id_k = m
	                        end
	                    end
	                    LightGraphs.add_edge!(g, l_all+1, id_k)
	                    LightGraphs.add_edge!(g, id_k, id_j)
	                    if lea_j[k] == true
	                        LightGraphs.add_edge!(g, id_j, id_k)
	                    end
	                end
				end

            end

        end



    end



	LightGraphs.collect(LightGraphs.edges(g))
	A = LightGraphs.LinAlg.adjacency_matrix(g)


	s_A = size(A, 1)
	sinks = []


	for i = 1:s_A
		num = 0
		if !(i in A.rowval)
			sink = []
			if i > length(all_cms)
				continue
			end
			for j in all_cms[i]
            	push!(sink, c_p_s[j])
			end
			push!(sinks, sink)
		end
	end

	return sinks
end

function sink_stuck_check(c_ids, sinks)
	perfect_insertion = []
	for id in c_ids
		push!(perfect_insertion, [id.p_id, id.e_id])
	end

	stuck_sink = []
	for sink in sinks
		sink_id = []
		for s_id in sink
			push!(sink_id, [s_id.p_id, s_id.e_id])
		end

		if !issubset(sink_id, perfect_insertion)
			push!(stuck_sink, sink)
		end
	end
	return stuck_sink
end


function is_goal_cm(c_ids, cm)
	perfect_insertion = []
	for id in c_ids
		push!(perfect_insertion, [id.p_id, id.e_id])
	end

	cm_ids = []
	for cp in cm
		push!(cm_ids, [cp.p_id, cp.e_id])
	end
	# println("here:",cm)
	if issubset(cm_ids, perfect_insertion)
		return true
	else
		return false
	end

end

function desired_sinks(c_ids, sinks)
	perfect_insertion = []
	for id in c_ids
		push!(perfect_insertion, [id.p_id, id.e_id])
	end

	desired_sinks = []
	for sink in sinks
		sink_id = []
		for s_id in sink
			push!(sink_id, [s_id.p_id, s_id.e_id])
		end

		if issubset(sink_id, perfect_insertion)
			push!(desired_sinks, sink)
		end
	end
	return desired_sinks
end

function insertion_error(peg, socket, sinks::Array)
	configs = []
	if length(sinks) == 0
		# println("zero length sinks")
		return configs
	end
	for sink in sinks
		cm = myType.cm_id(sink)
		result, c_list= ContactMode.is_cm_valid(peg, socket, cm)
		
		if result
			sort!(c_list, by = x->abs(x.th), rev=true)
			push!(configs, c_list[1].th)
		else
			# println("sinks are not valid")
			
		end
	end
	return configs
end

function sink_with_maximum_error(peg, socket, sinks::Array)
	configs = []
	for sink in sinks
		cm = myType.cm_id(sink)
		result, c_list= ContactMode.is_cm_valid(peg, socket, cm)
		if result
			sort!(c_list, by = x->abs(x.th), rev=true)
			push!(configs, abs(c_list[1].th))
		else
			# println("sinks are not valid")
			return nothing
		end
	end
	if length(configs) == 0
		return nothing
	end
	a = findmax(configs)
	return sinks[a[2]]
end

function sink_with_maximum_error_output_error(peg, socket, sinks::Array)
	configs = []
	for sink in sinks
		cm = myType.cm_id(sink)
		result, c_list= ContactMode.is_cm_valid(peg, socket, cm)
		if result
			sort!(c_list, by = x->abs(x.th), rev=true)
			push!(configs, abs(c_list[1].th))
		else
			# println("sinks are not valid")
			return nothing
		end
	end
	a = findmax(configs)
	# println("configs:", configs)
	return a[1]
end

function insertion_error(peg, socket, cm::myType.cm_id)
	config = 0

	result, c_list= ContactMode.is_cm_valid(peg, socket, cm)
	if result
		sort!(c_list, by = x->abs(x.th), rev=true)
		config = c_list[1].th
	else
		# println("cm is not valid: ")
		# println(cm)
		Draw.draw_sth(peg)
		Draw.draw_open_sth(socket)
		return nothing
	end

	return config
end


function test(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    # cp_ids, cms, c_lists = ContactMode.all_combinations(peg,socket,c_ids)

    # println(cms)
    sinks = get_sinks(peg,socket,c_ids)
	stuck_sink = sink_stuck_check(c_ids, sinks)
	configs = insertion_error(peg, socket, stuck_sink)
	# println(stuck_sink)
	# println(configs)
end

function test3(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    # cp_ids, cms, c_lists = ContactMode.all_combinations(peg,socket,c_ids)

    # println(cms)
	t = Adjust.compute_t(peg, socket)
	new_socket = Rotation.adj_upsocket(socket, 0.6)
	new_peg = Adjust.get_peg_from_t(peg, new_socket, t)
    cmt_graph(new_peg,new_socket,c_ids)

end

function test2(i)
    peg, socket = Construct_Peg_Socket.read_from_file(i)
    c_ids = Construct_Peg_Socket.get_contact_ids(i)
    es = mySocketEdges(socket)
    cp1 = contact_pair(1, peg.points[1], es.edges[1])
    cp2 = contact_pair(2, peg.points[2], es.edges[2])
	cp3 = contact_pair(2, peg.points[3], es.edges[2])
    cm1 = contact_mode([cp1])
    cm2 = contact_mode([cp1, cp2])
	cm3 = contact_mode([cp1, cp2, cp3])
	r1 = two_p_one_e(cm3)
	r2 = two_p_one_e(cm2)
	println(r1, r2)
    # CM_tansition_with_one_CP(cm1, cm2)
end


function test_initial_state(peg, socket, config)
    x_left = socket.points[1].x
    x_right = socket.points[length(socket.points)].x
    new_peg = Kinematics.get_peg_from_config(peg, config)
    p_x = new_peg.points[Int(floor((length(new_peg.points)+1)/2))].x
    # println("peg: ", new_peg)
    println("values: ", px, ", ", x_left, ", ", x_right, ", ", config.th)
    if p_x >= x_left && p_x <= x_right 
    	# need further testing, to make sure
    	# this would have a positive projection downwards

    	ang = Kinematics.get_max_socket_angle(socket)
    	println("max socket edge: ", ang)
    	t_a = pi/2 - (ang + config.th)
    	if tan(t_a) > CONST.friction_coeff
    		return true
    	end
    	return false
        # return true
    # elseif  p_x < x_left && config.th > CONST.friction_coeff/2
    #     return true
    # elseif p_x > x_right && config.th < CONST.friction_coeff/2
    #     return true
    else
        return false
    end
end


end






