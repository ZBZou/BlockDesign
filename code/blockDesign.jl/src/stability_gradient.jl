module Stability_gradient
import LinearAlgebra
import JuMP

import MathOptInterface
import NLopt
import Statistics


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket
import ..Insertion
import ..Convexhull
import ..Partial_order
import ..Rotation
import ..Adjust
import ..Insertion_gradient
import ..CMT_graph
import ..Draw

import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots




function all_close(a, b)
	if abs(a - b) < 0.011
		return true
	end
	return false;

end

function compute_signed_angle(v1, v2)

	vn = [0,0,1]


	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return th


end

function sink_stability(peg, socket, sink, dth = 0.01)
	f_cone = Insertion_gradient.friction_cone(socket)

	th_min = -pi
	th_max = pi

	while sink_can_move_under_th(peg, socket, f_cone, sink, th_min)
		th_min += dth
	end

	while sink_can_move_under_th(peg, socket, f_cone, sink, th_max)
		# println(th_max)
		th_max -= dth
	end

	# println("here: ", [th_min, th_max])
	return [th_min, th_max]
end

function sink_can_move_under_th(peg, socket, f_cone, sink, th)
	Fs = Insertion_gradient.force_analysis_for_sink(socket, peg, sink, th)
	n = length(f_cone)
	for i =1:length(sink)
		for j = 1:n
			if j == sink[i].e_id
				# println(f_cone[j], Fs[i])
				if Fs[i]>=f_cone[j][1] && Fs[i]<=f_cone[j][2]
					return false
				end
			end
		end
	end

	return true
end

function sinks(args)
	body
end

function sink_direction(peg, socket, sink)
	n = length(socket.points) - 1
	e_ids = []

	for cp in sink
		push!(e_ids, cp.e_id)
	end
	sort!(e_ids)

	if abs(e_ids[1] - 1) < abs(e_ids[length(e_ids)]-n)
		return 1
	elseif abs(e_ids[1] - 1) > abs(e_ids[length(e_ids)]-n)
		return -1
	else
		return 0
	end

end

function total_stab(peg, socket, sinks)
	n = length(sinks)
	left_stab_min = []
	left_stab_max = []
	right_stab_min = []
	right_stab_max = []
	for i = 1:n
		stab_i = sink_stability(peg, socket, sinks[i])
		r = sink_direction(peg, socket, sinks[i])
		if r == 1
			push!(right_stab_min, stab_i[1])
			push!(right_stab_max, stab_i[2])
		elseif r == -1
			push!(left_stab_min, stab_i[1])
			push!(left_stab_max, stab_i[2])
		else
			push!(right_stab_min, stab_i[1])
			push!(right_stab_max, stab_i[2])
			push!(left_stab_min, stab_i[1])
			push!(left_stab_max, stab_i[2])
		end

	end
	sort!(left_stab_min)
	sort!(left_stab_max, rev=true)
	sort!(right_stab_min)
	sort!(right_stab_max, rev=true)
	if length(left_stab_min) == 0
		return [0, 0]
	end

	if abs((abs(left_stab_min[1] - left_stab_max[1])) -  (abs(right_stab_min[1] - right_stab_max[1]))) > 0.01
		# println("is not symmetric")
		return [0,0]
	end

	return [left_stab_min[1], left_stab_max[1]]
end

function sink_with_least_stab(peg::myType.myPeg, socket::myType.mySocket,sinks::Array)
	n = length(sinks)
	stab = 100
	sink = []

	for i = 1:n
		stab_i = sink_stability(peg, socket, sinks[i])
		stab_value = abs(stab_i[1] - stab_i[2])
		if i == 1
			stab = stab_value
			sink = sinks[i]
			continue
		end
		if stab > stab_value
			stab = stab_value
			sink = sinks[i]
		end
	end

	return stab, sink
end



function stab_gradient(peg, socket, c_ids, ori_sinks, dth=0.1)

	n_s = length(socket.points)-1
	n_r = Int(floor(n_s/2))
	k_id = []

	r_mode = get_all_rotation_mode(socket)

	if r_mode == nothing
		# println("here")
		return nothing, nothing, nothing
	end

	# sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	#
	# stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
	# if stuck_sink != []
	# 	println("have stuck sinks")
	# 	return nothing, nothing, nothing
	# end
	#
	# desired_sinks = CMT_graph.desired_sinks(c_ids,sinks)
	# # stab = Statistics.mean(total_stab(peg, socket, desired_sinks))
	# stab, sink = sink_with_least_stab(peg, socket, desired_sinks)
	#
	# println("sink_with_least_stab: ", sink_with_least_stab)
	# # println("r_mode:", r_mode)

	# println(ori_sink)

	# stab_temp = sink_stability(peg, socket, ori_sink)
	
	if length(ori_sinks) == 0
		return nothing, nothing, 0
	end
	
	stab_temp = total_stab(peg, socket, ori_sinks)
	stab = abs(stab_temp[1] - stab_temp[2])

	new_peg = deepcopy(peg)
	new_socket = deepcopy(socket)


	for r in r_mode
		limit_r = edge_rotate_limit(socket, r)
		# println("limit_r:", limit_r)
		for i = 1:n_r

			if limit_r[i][1] || i == 1
				continue
			end

			k1 = 1
			k2 = -1
			socket1 = Rotation.rotate_socket(new_socket, r, i, dth*k1)
			peg1 = Adjust.get_adj_peg(new_peg, socket1, 0.)

			socket2 = Rotation.rotate_socket(new_socket, r, i, dth*k2)
			peg2 = Adjust.get_adj_peg(new_peg, socket2, 0.)

			# sinks1 = CMT_graph.get_sinks(peg1,socket1,c_ids)
			# stuck_sink1 = CMT_graph.sink_stuck_check(c_ids, sinks1)
			#
			# sinks2 = CMT_graph.get_sinks(peg2,socket2,c_ids)
			# stuck_sink2 = CMT_graph.sink_stuck_check(c_ids, sinks2)

			# stab1_temp = sink_stability(peg1, socket1, ori_sink)
			# stab1 = abs(stab1_temp[1] - stab1_temp[2])
			#
			# stab2_temp = sink_stability(peg2, socket2, ori_sink)
			# stab2 = abs(stab2_temp[1] - stab2_temp[2])
			stab1_temp = total_stab(peg1, socket1, ori_sinks)
			stab1 = abs(stab1_temp[1] - stab1_temp[2])

			stab2_temp = total_stab(peg2, socket2, ori_sinks)
			stab2 = abs(stab2_temp[1] - stab2_temp[2])

			# if stuck_sink1 == []
			# 	# stab1 = Statistics.mean(total_stab(peg1, socket1, sinks1))
			# 	stab1_temp = sink_stability(peg1, socket1, sink)
			# 	stab1 = abs(stab1_temp[1] - stab1_temp[2])
			# else
			# 	stab1 = 0
			# end
			#
			# if stuck_sink2 == []
			# 	# stab2 = Statistics.mean(total_stab(peg2, socket2, sinks2))
			# 	stab2_temp = sink_stability(peg2, socket2, sink)
			# 	stab2 = abs(stab2_temp[1] - stab2_temp[2])
			# else
			# 	stab2 = 0
			# end

			# PyPlot.figure()
			# Draw.draw_open_sth(new_socket)
			# Draw.draw_open_sth(socket1, (0,0,1))
			#
			# PyPlot.figure()
			# Draw.draw_open_sth(new_socket)
			# Draw.draw_open_sth(socket2, (0,1,1))

			# println([stab, stab1, stab2])
			if stab >= stab1 && stab >= stab2
				push!(k_id, [stab, r, i, 0])
				#return nothing, nothing
			elseif stab1 > stab2
				push!(k_id, [stab1, r, i, k1])
				#return peg1, socket1
			else
				push!(k_id, [stab2, r, i, k2])
				#return peg2, socket2
			end
		end
	end

	if k_id != []

		sort!(k_id, by = x->x[1], rev = true)
		# println("k_id:", k_id)

		r = k_id[1][2]
		index = k_id[1][3]
		k = k_id[1][4]

		return r, index, k
	end

	return nothing, nothing, 0

end

function is_sinks_equal(new_sinks,sinks)
	n1 = length(new_sinks)
	n2 = length(sinks)

	if n1 != n2
		return false
	end

	for i = 1:n1
		for j = 1:length(sinks[i])
			if new_sinks[i][j].p_id != sinks[i][j].p_id || new_sinks[i][j].e_id != sinks[i][j].e_id
				return false
			end
		end
	end

	return true
end


function socket_optimize(peg, socket, c_ids, i = [], pic_count = 0, dth=0.1, k_dt = 0.01)
	sinks = CMT_graph.get_sinks(peg,socket,c_ids)
	stuck_sink = CMT_graph.sink_stuck_check(c_ids, sinks)
	desired_sinks = CMT_graph.desired_sinks(c_ids, sinks)
	up_limit = Insertion_gradient.rotation_upper_limit(peg, socket, desired_sinks)
	e_angle = Insertion_gradient.entry_angle(socket)
	new_peg = deepcopy(peg)
	output_peg = deepcopy(peg)
	new_socket = deepcopy(socket)
	new_sinks = deepcopy(sinks)
	output_socket = deepcopy(socket)

	# ori_stab = total_stab(peg, socket, desired_sinks)
	# ori_stab, ori_sink = sink_with_least_stab(peg, socket, desired_sinks)

	r, index, k = stab_gradient(new_peg, new_socket, c_ids, desired_sinks, dth)



	while stuck_sink == [] && abs(e_angle[1]-e_angle[2]) < abs(up_limit[1]-up_limit[2]) && k!= 0 && is_sinks_equal(new_sinks,sinks)
		t = Adjust.compute_t(new_peg, new_socket)
		if k == nothing
			return peg, socket, pic_count
		end
		new_socket = Rotation.rotate_socket(new_socket, r, index, k_dt*k)
		new_peg = Adjust.get_adj_peg(new_peg, new_socket, 0.)
		r, index, k = stab_gradient(new_peg, new_socket, c_ids, new_sinks, dth)
		new_sinks = CMT_graph.get_sinks(new_peg,new_socket,c_ids)
		# println("here1,k: ", k)
		stuck_sink = CMT_graph.sink_stuck_check(c_ids, new_sinks)
		# println("here2 ")
		desired_sinks = CMT_graph.desired_sinks(c_ids, new_sinks)
		# println("here3 ")
		up_limit = Insertion_gradient.rotation_upper_limit(new_peg, new_socket, desired_sinks)
		# println("here4 ")
		e_angle = Insertion_gradient.entry_angle(new_socket)

		if stuck_sink != [] && !is_sinks_equal(new_sinks,sinks)
			# println("here5 ")
			break
		end
		output_peg = deepcopy(new_peg)
		output_socket = deepcopy(new_socket)

		Draw.draw_anime(new_peg, new_socket, i, pic_count)
		# PyPlot.clf()
		# PyPlot.axis("off")
		# Draw.draw_sth(output_peg, (1,0,0), 3.)
		# Draw.draw_open_sth(output_socket, (0.5,.4,.2), 3.)
		pic_count += 1
	end


	return output_peg, output_socket, pic_count

end

function peg_gradient(peg::myType.myPeg, socket::myType.mySocket,
	F=[0,0,0], closed = false, convex = false, limit_idex = [])
	#find the gradiant for peg design to stability, e.g., dt
	closed = true
	idset = []
	thset, fth, sth= Rotation.get_rotate_th(peg, socket, F, closed, convex)
	thsum = fth[1] + sth[1]
	# println("thsum:", thsum)
	#println("thset:", thset)
	idset = push!(idset, fth[2])
	idset = push!(idset, sth[2])

	n_t = Int(floor(length(peg.points)/2))
	#println(n_t)
	# if dt == []
	# 	for i = 1:n_t
	# 		push!(dt, 0.01)
	# 	end
	# end

	#linear optimization, thus test limit value
	dt = []
	for i = 1:n_t
		push!(dt, 1.)
	end



	g = zeros(n_t)
	k_new = zeros(n_t)
	k_id = Array{Array}(undef, 0)

	limit_t = Adjust.check_t_limit(peg, socket)
	for i = 1:n_t

		if limit_t[i] || i in limit_idex
			# println("limit_t: ", i)
			continue
		end

		g1 = deepcopy(g)
		g2 = deepcopy(g)
		g1[i] = 1
		g2[i] = -1

		peg1 = Adjust.get_adj_peg(peg, socket, g1.*dt)
		peg2 = Adjust.get_adj_peg(peg, socket, g2.*dt)

		thset1, fth1, sth1= Rotation.get_rotate_th(peg1, socket, F, closed, convex)
		thset2, fth2, sth2= Rotation.get_rotate_th(peg2, socket, F, closed, convex)
		# thsum1 = 0
		# count1 = 0
		# thsum2 = 0
		# count2 = 0
		# for th1 in thset1
		# 	if th1[2] in idset
		# 		thsum1 += th1[1]
		# 		count1 += 1
		# 	end
		# end
		# if count1 != length(idset)
		# 	thsum1 = fth1[1] + sth1[1]
		# end
		#
		#
		#
		# for th2 in thset2
		# 	if th2[2] in idset
		# 		thsum2 += th2[1]
		# 		count2 += 1
		# 	end
		# end
		#
		# if count2 != length(idset)
		# 	thsum2 = fth2[1] + sth2[1]
		# end

		# if MoC change, CMT must change, still compare max rotation angle
		thsum1 = 0
		thsum2 = 0

		thsum1 = fth1[1] + sth1[1]
		thsum2 = fth2[1] + sth2[1]

		# println("thset1:", thset1)
		# println("thset2:", thset2)

		 # println("thsum1:", thsum1)
		 # println("thsum2:", thsum2)


		if thsum <= thsum1 && thsum <= thsum2

			push!(k_id, [thsum, i, 0])
		elseif thsum1 < thsum2

			push!(k_id, [thsum1, i, 1])
		else

			push!(k_id, [thsum2, i, -1])
		end
	end




	if k_id != []
		sort!(k_id, by = x->x[1])
		id = Int(k_id[1][2])
		k_new[id] = k_id[1][3]
	end

	# println("peg k_id : ", k_id)
	return k_new


end

function Permutation(N::Int, A = 2)

	return reverse.(digits.(0:A^N-1,base=A,pad=N))
end



function get_all_rotation_mode(socket::myType.mySocket)
	n_s = length(socket.points) - 1
	n_r = Int(floor(n_s/2))

	if n_s <= 1
		return nothing
	end


	r_mode = Permutation(n_r)

	if length(r_mode) == 1
		return r_mode
	end
	count = 1
	for r in r_mode
		if length(r) < 2
			return r_mode[count]
		end
		count += 1
		r[1] = 0
		r[2] = 0
		if iseven(n_s)

			if n_r == 2
				return nothing
			end
			r[length(r)] = 1
		end
	end

	r_mode = unique(r_mode)

	return r_mode

end

function test_all_rmode(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	r = get_all_rotation_mode(socket)
	if r != nothing
		# println(r)
	else
		show(r)
	end
end



function edge_rotate_limit(socket::myType.mySocket, r, ep = 0.005)
	n_s = length(socket.points) - 1
	n_r = length(r)
	limit = []

	for i = 1:n_r
		p_id = i
		if r[i] == 0
			p1 = socket.points[p_id]
			p2 = socket.points[p_id+1]
			p3 = socket.points[p_id+2]
			v1 = [p2.x - p1.x, p2.y - p1.y, 0]
			v2 = [p3.x - p1.x, p3.y - p1.y, 0]
			th1 = compute_signed_angle(v1, [1, 0, 0])
			th2 = compute_signed_angle(v2, [1, 0, 0])
			th_test = compute_signed_angle(v1, v2)

			# if th_test < 0, concave

			if  i == n_r && iseven(n_s)
				mid = [(p3.x+p2.x)/2, p2.y]
				temp = [mid[1] - p1.x, mid[2] - p1.y, 0]
				minth = compute_signed_angle(temp, [1, 0, 0])
				maxth = pi/2
			elseif th_test <= 0

				minth = 0
				maxth = th2
			else
				minth = th2
				maxth = pi/2
			end

			# println([th_test, th1, maxth])
			if th1 - ep <= minth || th1 + ep >= maxth
				push!(limit, [true, r, i])

			else
				push!(limit, [false, r, i])

			end
		else
			if i == 1 || i == 2
				error("fix socket first edge")
			end

			p1 = socket.points[p_id+1]
			p2 = socket.points[p_id]
			p3 = socket.points[p_id-1]
			v1 = [p2.x - p1.x, p2.y - p1.y, 0]
			v2 = [p3.x - p1.x, p3.y - p1.y, 0]
			th1 = compute_signed_angle(v1, [-1, 0, 0])
			th2 = compute_signed_angle(v2, [-1, 0, 0])
			th_test = compute_signed_angle(v1, v2)


			if th_test >= 0
				minth = th2
				maxth = pi/2
			else
				minth = 0
				maxth = th2
			end
			# println(th_test,",", th1)
			# println("minth:", minth)
			# println("maxth:", maxth)


			if th1 + ep >= maxth || th1 - ep <= minth
				push!(limit, [true, r, i])

			else
				push!(limit, [false, r, i])

			end
		end
	end

	# fix first edge
	# limit[1][1] = true

	return limit

end

function socket_limit(socket)



end


function socket_gradient(peg::myType.myPeg, socket::myType.mySocket,
	F=[0,0,0], convex = false, dth=0.005)
	Insertion_closed = false
	Stability_closed = true
	n_s = length(socket.points) - 1
	n_r = Int(floor(n_s/2))

	# 0 means rotate point is start point, 1 means rotate point is end point
	# first and last edges should be rotated by the start point
	# as the height of the socket should be unchanged, if n_s is odd, the mid two edges should only be rotated by end point
	r_mode = get_all_rotation_mode(socket)

	if r_mode == nothing
		return nothing, nothing, nothing
	end

	k_id = []
	idset = []
	thset, fth, sth= Rotation.get_rotate_th(peg, socket, F, Stability_closed, convex)
	thsum = fth[1] + sth[1]
	#println("thset:", thset)
	idset = push!(idset, fth[2])
	idset = push!(idset, sth[2])
	# println("thsum:", thsum)
	#th_limit = Insertion_gradient.limit_th_with_x_axis(peg::myType.myPeg, socket::myType.mySocket, F, Insertion_closed, convex)

	# println("r_mode:", r_mode)

	for r in r_mode
		limit_r = edge_rotate_limit(socket, r)
		# println("limit_r:", limit_r)
		for i = 1:n_r

			if limit_r[i][1] || i == 1
				continue
			end

			k1 = 1
			k2 = -1
			socket1 = Rotation.rotate_socket(socket, r, i, dth*k1)
			peg1 = Adjust.get_adj_peg(peg, socket1, 0.)

			socket2 = Rotation.rotate_socket(socket, r, i, dth*k2)
			peg2 = Adjust.get_adj_peg(peg, socket2, 0.)

			thset1, fth1, sth1= Rotation.get_rotate_th(peg1, socket1, F, Stability_closed, convex)
			thset2, fth2, sth2= Rotation.get_rotate_th(peg2, socket2, F, Stability_closed, convex)
			if thset1 == nothing || thset2 == nothing
				return nothing, nothing, nothing
			end



			#println("thset2:", thset2)
			# thsum1 = 0
			# count1 = 0
			# thsum2 = 0
			# count2 = 0
			# for th1 in thset1
			# 	if th1[2] in idset
			# 		thsum1 += th1[1]
			# 		count1 += 1
			# 	end
			# end
			# if count1 != length(idset)
			# 	thsum1 = fth1[1] + sth1[1]
			# end
			#
			#
			#
			# for th2 in thset2
			# 	if th2[2] in idset
			# 		thsum2 += th2[1]
			# 		count2 += 1
			# 	end
			# end
			#
			# if count2 != length(idset)
			# 	thsum2 = fth2[1] + sth2[1]
			# end
			# thsum1 = fth1[1] + sth1[1]
			# thsum2 = fth2[1] + sth2[1]

			thsum1 = 0
			thsum2 = 0

			thsum1 = fth1[1] + sth1[1]
			thsum2 = fth2[1] + sth2[1]

			 # println("thsum1:",thsum1)
			 # println("thsum2:",thsum2)


			if thsum <= thsum1 && thsum <= thsum2
				push!(k_id, [thsum, r, i, 0])
				#return nothing, nothing
			elseif thsum1 < thsum2
				push!(k_id, [thsum1, r, i, k1])
				#return peg1, socket1
			else
				push!(k_id, [thsum2, r, i, k2])
				#return peg2, socket2
			end
		end
	end

	if k_id != []

		sort!(k_id, by = x->x[1])
		# println("k_id:", k_id)

		r = k_id[1][2]
		index = k_id[1][3]
		k = k_id[1][4]

		return r, index, k
	end

	return nothing, nothing, nothing

end

function test_socketgradient(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	socket_gradient(peg, socket)


end

function test_ifchangeth(i)
	dth = 0.005
	peg, socket = Construct_Peg_Socket.read_from_file(i)

	r, index, k = Stability_gradient.socket_gradient(peg, socket)
	#println("r, index, k : ", r, index, k)

	 socket_new = Rotation.rotate_socket(socket, r, index, k*dth)
	 peg_new = Adjust.get_adj_peg(peg, socket_new, 0.)

	r, index, k = Stability_gradient.socket_gradient(peg_new, socket_new)

	#thset1, fth1, sth1= Rotation.get_rotate_th(peg1, socket1, F, Stability_closed, convex)





end

function socket_gradient_low(peg::myType.myPeg, socket::myType.mySocket,
	F=[0,0,0], closed = false, convex = false, dth=0.005)
	#find the gradiant for low edges of socket to stability,
	#e.g., the edge rotation direction
	closed = true

	idset = []
	thset, fth, sth= Rotation.get_rotate_th(peg, socket, F, closed, convex)
	thsum = fth[1] + sth[1]
	# println(fth,",", sth)


	idset = push!(idset, fth[2])
	idset = push!(idset, sth[2])
	k1 = 1
	k2 = -1

	socket1 = Rotation.adj_lowsocket(socket, dth*k1)
	peg1, full1 = Adjust.get_adj_peg(peg, socket1, 0.)

	socket2 = Rotation.adj_lowsocket(socket, dth*k2)
	peg2, full2 = Adjust.get_adj_peg(peg, socket2, 0.)


	thset1, fth1, sth1= Rotation.get_rotate_th(peg1, socket1, F, closed, convex)
	thset2, fth2, sth2= Rotation.get_rotate_th(peg2, socket2, F, closed, convex)
	thsum1 = 0
	thsum2 = 0
	for th1 in thset1
		if th1[2] in idset
			thsum1 += th1[1]
		end
	end

	for th2 in thset2
		if th2[2] in idset
			thsum2 += th2[1]
		end
	end

	# println(fth1,",", sth1)
 	# println(fth2,",", sth2)


	if thsum1 >= thsum && thsum2 >= thsum
		return 0
	elseif thsum1 <= thsum2
		return k1
	else
		return k2

	end


end


function socket_gradient_up(peg::myType.myPeg, socket::myType.mySocket,
	F=[0,0,0], closed = false, convex = false, dth=0.005)
	#find the gradiant for up edges of socket to stability,
	#e.g., the edge rotation direction
	closed = true

	idset = []
	thset, fth, sth= Rotation.get_rotate_th(peg, socket, F, closed, convex)
	thsum = fth[1] + sth[1]
	# println(fth,",", sth)


	idset = push!(idset, fth[2])
	idset = push!(idset, sth[2])
	k1 = 1
	k2 = -1

	socket1 = Rotation.adj_upsocket(socket, dth*k1)
	peg1, full1 = Adjust.get_adj_peg(peg, socket1, 0.)

	socket2 = Rotation.adj_upsocket(socket, dth*k2)
	peg2, full2 = Adjust.get_adj_peg(peg, socket2, 0.)


	thset1, fth1, sth1= Rotation.get_rotate_th(peg1, socket1, F, closed, convex)
	thset2, fth2, sth2= Rotation.get_rotate_th(peg2, socket2, F, closed, convex)
	thsum1 = 0
	thsum2 = 0
	for th1 in thset1
		if th1[2] in idset
			thsum1 += th1[1]
		end
	end

	for th2 in thset2
		if th2[2] in idset
			thsum2 += th2[1]
		end
	end

	# println(fth1,",", sth1)
 	# println(fth2,",", sth2)


	if thsum1 >= thsum && thsum2 >= thsum
		return 0
	elseif thsum1 <= thsum2
		return k1
	else
		return k2

	end


end


function test(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
    F = [0,0,0]
    closed = true
    convex = Insertion.is_convex(socket)

	peg_gradient(peg, socket, F, closed, convex)

end

function test_low(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
    F = [0,0,0]
    closed = true
    convex = Insertion.is_convex(socket)

	socket_gradient_low(peg, socket, F, convex)

end

function test_up(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
    F = [0,0,0]
    closed = true
    convex = Insertion.is_convex(socket)

	socket_gradient_up(peg, socket, F, convex)

end















end
