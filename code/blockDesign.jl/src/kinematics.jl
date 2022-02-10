module Kinematics # forward kinematics

import ..myType

function rotation_matrix_2D(angle::Float64)
	R = Matrix{Float64}(undef, 2, 2)

	R[1, 1] = cos(angle)
	R[1, 2] = - sin(angle)
	R[2, 1] = sin(angle)
	R[2, 2] = cos(angle)

	return R
end

function transformation_matrix_2D(x::Float64, y::Float64, theta::Float64)
	R = Matrix{Float64}(undef, 3, 3)

	R[1, 1] = cos(theta)
	R[1, 2] = - sin(theta)
	R[1, 3] = x
	R[2, 1] = sin(theta)
	R[2, 2] = cos(theta)
	R[2, 3] = y
	R[3, 1] = 0;
	R[3, 2] = 0
	R[3, 3] = 1

	return R
end

function transform_indexed_point(R::Matrix{Float64}, point::myType.Indexed_point)

	v = Array{Float64}(undef, 3)
	v[1] = point.x
	v[2] = point.y
	v[3] = 1

	result = R * v

	return myType.Indexed_point(point.index, result[1], result[2])

end

function transform_point(R::Matrix{Float64}, p_x, p_y)

	v = Array{Float64}(undef, 3)
	v[1] = p_x
	v[2] = p_y
	v[3] = 1

	result = R * v

	return [result[1], result[2]]

end

function compute_points_from_config(peg::myType.myPeg, config::myType.configSE2)
	# myPeg points are in peg frame

	ps = Array{myType.Indexed_point}(undef, size(peg.points, 1))

	R = transformation_matrix_2D(config.x, config.y, config.th)

	for i = 1:size(peg.points, 1)

		ps[i] = transform_indexed_point(R, peg.points[i])

	end

	return ps

end

function get_peg_from_config(peg::myType.myPeg, config::myType.configSE2)
	ps = compute_points_from_config(peg, config)


	return myType.myPeg(ps)
end

function get_pegs_from_configs(peg::myType.myPeg, configs::Array)
	n = length(configs)
	pegs = []
	for i = 1:n
		push!(pegs, get_peg_from_config(peg, configs[i][1]))
    end


	return pegs
end

function rotate_peg_around_p(peg::myType.myPeg, p2::myType.myPoint2d, theta::Float64)
	n = length(peg.points)

	for i = 1:n
		peg.points[i] = rotate_p1_around_p2(peg.points[i], p2, theta)
	end


	return peg
end

function rotate_peg_around_p(peg::myType.myPeg, p2::myType.Indexed_point, theta::Float64)
	n = length(peg.points)

	for i = 1:n
		peg.points[i] = rotate_p1_around_p2(peg.points[i], p2, theta)
	end


	return peg
end

function rotate_p1_around_p2(p1::myType.Indexed_point, p2::myType.myPoint2d, theta::Float64)



	R = rotation_matrix_2D(theta)

	v1 = [p1.x, p1.y]
	v2 = [p2.x, p2.y]
	v = v2+R*(v1-v2)
	return myType.Indexed_point(p1.index, v[1], v[2])
	# x = cos(theta) * (p1.x-p2.x) - sin(theta) * (p1.y-p2.y) + p2.x
	#
	# y = sin(theta) * (p1.x-p2.x) + cos(theta) * (p1.y-p2.y) + p2.y
	# println(x, y)
	#
	# return myType.Indexed_point(p1.index, x, y)
end

function rotate_p1_around_p2(p1::myType.Indexed_point, p2::myType.Indexed_point, theta::Float64)
	x = cos(theta) * (p1.x-p2.x) - sin(theta) * (p1.y-p2.y) + p2.x

	y = sin(theta) * (p1.x-p2.x) + cos(theta) * (p1.y-p2.y) + p2.y

	return myType.Indexed_point(p1.index, x, y)
end

function rotate_p1_around_p2(p1::Array, p2::Array, theta::Float64)
	x = cos(theta) * (p1[1]-p2[1]) - sin(theta) * (p1[2]-p2[2]) + p2[1]

	y = sin(theta) * (p1[1]-p2[1]) + cos(theta) * (p1[2]-p2[2]) + p2[2]

	return [x, y]
end

# function get_config_from_mode(peg::myPeg, c_m::contact_mode)

# end

# function get_config_constraint_from_pair(peg::myPeg, c_p::contact_pair)



# end


function get_max_socket_angle(socket) 
	max_ang = 0
	len = Int(floor(length(socket.edges)/2))
	for i = 1: len
		e = socket.edges[i]
		p_i = socket.points[e.s]
		p_j = socket.points[e.e]
		ang = abs(atan((p_j.x - p_i.x) / (p_j.y - p_i.y)))
		# println("pi: ", p_i, ", pj: ", p_j)
		# println("ang: ", atan((p_j.x - p_i.x) / (p_j.y - p_i.y)))
		if ang > max_ang


			max_ang = ang
		end
	end
	return max_ang
end







export compute_points_from_config, rotation_matrix_2D, get_peg_from_config

end
