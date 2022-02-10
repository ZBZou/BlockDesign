module Adjust
import LinearAlgebra
import JuMP

import MathOptInterface
import NLopt


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket
import ..Draw



import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots

function all_close(a, b)
	if abs(a - b) < 0.02
		return true
	end
	return false;

end

function get_adj_peg(peg::myType.myPeg, socket::myType.mySocket, dt=0.1, limit = 0.95)

	t = compute_t(peg, socket)
	te = ep_trans_t(socket, t)
	# println(t, ",", dt)

	for i = 1:length(t)
		t[i] += dt
	end


	t = bound_t(t, te, limit)

	new_peg = get_peg_from_t(peg, socket, t)

	return new_peg


end

function bound_t(t::Array{Float64}, te::Array, limit::Float64)

	for i = 1 : length(t)
		if length(te[i]) == 1
			if t[i] <= te[i]+(1-limit)
				t[i] = te[i]+(1-limit)
			elseif t[i] >= limit-te[i]
				t[i] = limit-te[i]
			end
		else
			if t[i] <= te[i][1]+(1-limit)
				t[i] = te[i][1]+(1-limit)
			elseif t[i] >= limit-te[i][2]
				t[i] = limit-te[i][2]
			end
		end

	end
	# println("update_t:", t)
	return t

end

function check_t_limit(peg::myType.myPeg, socket::myType.mySocket, limit = 0.95)
	result = []
	t = compute_t(peg, socket)
	te = ep_trans_t(socket, t)

	n_p = length(peg.points)
	n_s = length(socket.points) - 1


	for i = 1 : length(t)
		if length(te[i]) == 1
			if t[i] <= te[i]+(1-limit)
				t[i] = te[i]+(1-limit)
				push!(result, true)
			elseif t[i] >= limit+2*te[i]
				t[i] = limit+2*te[i]
				push!(result, true)
			else
				push!(result, false)
			end
		else
			if t[i] <= te[i][1]+(1-limit)
				t[i] = te[i][1]+(1-limit)
				push!(result, true)
			elseif t[i] >= limit-te[i][2]
				t[i] = limit-te[i][2]
				push!(result, true)
			else
				push!(result, false)
			end
		end


	end
	# println(result)
	return result
end

function check_peg_limit(peg::myType.myPeg, socket::myType.mySocket, limit = 0.95)
	result = check_t_limit(peg, socket, limit)
	peg_limit = true

	for i in result
		if !i
			peg_limit = false
		end
	end

	return peg_limit



end


function get_adj_peg(peg::myType.myPeg, socket::myType.mySocket, dt::Array{Float64}, limit = 0.99)

	t = compute_t(peg, socket)
	te = ep_trans_t(socket, t)

	# Draw.draw_sth(peg)
	# Draw.draw_open_sth(socket)
	# println(t, ",", dt)
	t = t+dt

	t = bound_t(t, te, limit)

	new_peg = get_peg_from_t(peg, socket, t)

	return new_peg


end

function test_pegadjust(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	n = Int(floor(length(peg.points)/2))
	dt = Array{Float64}(undef, 0)
	for i = 1:n
		push!(dt, 0.)
	end

	new_peg = get_adj_peg(peg, socket, 0.2)

	Draw.draw_open_sth(socket)
	Draw.draw_sth(peg)
	Draw.draw_sth(new_peg)

end

function compute_t(peg::myType.myPeg, socket::myType.mySocket)
	n_s = length(socket.edges) - 1
	n_p = length(peg.points)
	t = Array{Float64}(undef, 0)
	# if n_s > 6 || n_s < 3 || abs(n_s - n_p) >=2
	# 	error("not ideal joint")
	# end

	n = Int(floor(n_p/2))

	for i = 1:n
		p_i = socket.points[i]
		p_j = socket.points[i+1]
		p_p = peg.points[i]
		v1 = [p_j.x - p_i.x, p_j.y - p_i.y]
		v2 = [p_p.x - p_i.x, p_p.y - p_i.y]
		ti = v_ratio(v1, v2)
		push!(t, ti)
	end
	
	# if n_p % 2 == 1
	# 	push!(t, 0.5)
	# end


	return t

end

function test_t(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	t = compute_t(peg, socket)


end

function ep_trans_t(socket::myType.mySocket, t::Array, ep=0.1)
	te = []
	n_s = length(socket.edges)

	n = length(t)

	for i = 1:n

		p_i = socket.points[i]
		p_j = socket.points[i+1]
		v = [p_j.x - p_i.x, p_j.y - p_i.y, 0]
		l = vector_length(v)

		th = compute_abs_angle(v, [1,0,0])

		# if n == n_s/2  && iseven(n_s) && i == n
		if i == n
			lep1 = ep*tan(pi/2 - th)
			t1 = lep1/l
			lep2 = ep*cos(pi/2 - th)
			t2 = lep2/l
			push!(te, [t1, t2])

		else

			lep = ep*tan(pi/2 - th)
			tei = lep/l
			push!(te, tei)
		end
	end

	return te
end

function test_te(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	te = ep_trans_t(socket, [0.4,0.4])

end




function get_peg_from_t(peg::myType.myPeg, socket::myType.mySocket, t::Array{Float64},
	 ep=0.1)
	
	ep = 0.05
	 
	n_p = length(peg.points)
	n_s = length(socket.edges)-1
	# assume ideal peg are symmetric.

	mid_x = (socket.points[1].x + socket.points[length(socket.points)].x) /2

	new_peg = deepcopy(peg)
	# if n_s > 6 || n_s < 3 || abs(n_s - n_p) >=2
	# 	error("not ideal joint")
	# end
	# println("length: ", length(t))

	for i = 1 : length(t)
		p_i = socket.points[i]
		p_j = socket.points[i+1]
		p = compute_xy(p_i, p_j, t[i], ep)
		new_peg.points[i].x = p[1]
		new_peg.points[i].y = p[2]
		new_peg.points[n_p+1-i].x = mid_x-p[1]
		new_peg.points[n_p+1-i].y = p[2]
	end

	


	return new_peg

end


function compute_signed_angle(v1, v2)

	vn = [0,0,1]

	#os = LinearAlgebra.dot(v1,v2)/v1[1]^2+v2

	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return th


end

function v_ratio(v1, v2)
	l = vector_length(v1)
	return abs((v1[1]*v2[1]+v1[2]*v2[2])/l^2)

end


function vector_length(v::Array)
	scale = 0
	for i = 1:length(v)
		scale += v[i]*v[i]
	end

	return sqrt(scale)

end


function compute_abs_angle(v1, v2)

	vn = [0,0,1]

	#os = LinearAlgebra.dot(v1,v2)/v1[1]^2+v2

	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return abs(th)


end


function compute_xy(p1::myType.Indexed_point, p2::myType.Indexed_point, t::Float64,
	ep::Float64)


	# v = [p2.x - p1.x, p2.y - p1.y, 0]
	# l = vector_length(v)
	# p_x = t*l
	# p_y = ep
	# th = compute_signed_angle([1,0,0], v)
	# R  = Kinematics.transformation_matrix_2D(p1.x, p1.y, th)
	# p_new = Kinematics.transform_point(R, p_x, p_y)
	# return p_new
	# 
	v = [p2.x-p1.x, p2.y-p1.y]
	
	n = [-v[2], v[1]]
	
	np = [p1.x + t * v[1], p1.y + t * v[2]]
	
	p_new = [np[1] + ep * n[1], np[2] + ep * n[2]]
	return p_new

end












end
