module CMT

import LinearAlgebra
import JuMP

import MathOptInterface
import NLopt


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket


#import ..myType
using ..myType

function angle_diff(ang1, ang2)

	# first convert all to 0 to 2pi
	while ang1 < 0
		ang1 += 2*pi
	end
	while ang2 < 0
		ang2 += 2*pi
	end

	while ang1 > 2*pi
		ang1 -= 2*pi
	end

	while ang2 > 2*pi
		ang2 -= 2 * pi
	end

	# now, compute the difference

	diff = abs(ang1 - ang2)
	if diff > pi
		diff = 2*pi - diff
	end

	return diff

end


function negate_vector(v::Array)
	for i = 1:length(v)
		v[i] = -v[i]
	end

end


function vector_length(v::Array)
	scale = 0
	for i = 1:length(v)
		scale += v[i]*v[i]
	end

	return sqrt(scale)

end



function check_violation_in_seq(peg::myPeg, config::configSE2, socket::mySocket,
			current::Int)

	left = false
	right = false



	ps = compute_points_from_config(peg, config)

	# the vio list needs to be an array of tuples
	# it contains the contact point and contact edge
	# or, can be just a contact_pair?
	# decision:
	# a pair of indecies would be sufficient
	# first, index of contact point, second, index of contact edge

	vio = Array{c_pair}(undef, 0)

	for i = 1:size(ps, 1)
		if i == current
			continue
		end

		result = Constraints.inside_socket(ps[i], socket)
		if result > 0
			# push!(vio, i)
			push!(vio, c_pair(i, result, 1))
			if i < current
				left = true
			elseif i > current
				right = true
			end
		end


	end




	# this only test if for a given configuration
	# whether the constraints are violated
	# need to find if exist or not a valid configuration

	return left, right, vio

end

function distance(p1::Indexed_point, p2::Indexed_point)
	return sqrt((p1.x-p2.x)^2 + (p1.y-p2.y)^2)
end

function distance(p1::Indexed_point, p2::myPoint2d)
	return sqrt((p1.x-p2.x)^2 + (p1.y-p2.y)^2)
end

function distance(p1::myPoint2d, p2::Indexed_point)
	return sqrt((p1.x-p2.x)^2 + (p1.y-p2.y)^2)
end


function get_line_eq(p1::myType.Indexed_point, p2::myType.Indexed_point)

	# use p1 and p2 to get a line equation in the form of
	# ax + by + c = 0


	# this way, the returned equation can be directly used in the constraint

	# set a = 1

	line = Array{Float64}(undef, 3)

	if all_close(p1.x, p2.x)
		line[1] = 1
		line[2] = 0
		line[3] = -p1.x
	elseif all_close(p1.y, p2.y)
		# second case, similar y, horizontal line
		line[1] = 0
		line[2] = 1
		line[3] = -p1.y
	else
		# third case, not special line
		# use x + by + c = 0
		line[1] = 1
		line[3] = (p2.x/p2.y - p1.x/p1.y) / (1/p1.y-1/p2.y)
		line[2] = (-p1.x-line[3]) / p1.y
	end

	return line

end

function line_circle_intersection(p1::myPoint2d, p2::myPoint2d, center::myPoint2d,
		radius::Float64)
	# current, assume that p1 to p2 is counter clockwise
	#

	# vector pointing to p2
	v = [p2.x-p1.x, p2.y-p1.y]

	# define t

	A = v[1]^2 + v[2]^2
	B = 2 * (v[1]*(p1.x-center.x)+v[2]*(p1.y-center.y))
	C = (p1.x-center.x)^2 + (p1.y-center.y)^2 - radius^2

	# above is the equation At^2 + Bt + C = 0;

	if B^2 - 4*A*C < 0
		# should not happen, but should test anyway for safety
		println("Error! line not intersecting with circle!")
		exit(1)
	end


	# find solution
	t1 = (-B + sqrt(B^2-4*A*C)) / (2 * A)
	t2 = (-B - sqrt(B^2-4*A*C)) / (2 * A)


	return t1, t2

end

function line_circle_intersection(p1::Indexed_point, p2::Indexed_point, center::myPoint2d,
		radius::Float64)
	# current, assume that p1 to p2 is counter clockwise
	#

	# vector pointing to p2
	v = [p2.x-p1.x, p2.y-p1.y]

	# define t

	A = v[1]^2 + v[2]^2
	B = 2 * (v[1]*(p1.x-center.x)+v[2]*(p1.y-center.y))
	C = (p1.x-center.x)^2 + (p1.y-center.y)^2 - radius^2

	# above is the equation At^2 + Bt + C = 0;

	if B^2 - 4*A*C < 0
		# should not happen, but should test anyway for safety
		println("Error! line not intersecting with circle!")
		exit(1)
	end


	# find solution
	t1 = (-B + sqrt(B^2-4*A*C)) / (2 * A)
	t2 = (-B - sqrt(B^2-4*A*C)) / (2 * A)


	return t1, t2

end


function angle_between_rotation(point::myPoint2d, before::myPoint2d, after::myPoint2d)
	v1 = [before.x - point.x, before.y - point.y]
	v2 = [after.x - point.x, after.y - point.y]

	nv1 = sqrt(v1[1]^2+v1[2]^2)
	nv2 = sqrt(v2[1]^2+v2[2]^2)

	return acos((LinearAlgebra.dot(v1, v2)) / (nv1*nv2))

end

function rotate_left_against_edge(peg::myPeg, socket::mySocket, current::Int, point::myPoint2d,
		contact_point::Int, contact_edge::Int)


	vio_point = peg.points[contact_point]
	vio_edge_id = socket.edges[contact_edge]

	vio_edge_pi = socket.points[vio_edge_id]
	vio_edge_pj = socket.points[vio_edge_id]

	dist = distance(point, vio_point)

	t1, t2 = line_circle_intersection(vio_edge_pi, vio_edge_pj, point, dist)

	v = [vio_edge_pj.x - vio_edge_pi.x, vio_edge_pj.y - vio_edge_pi.y]
	new_point = myPoint2d(vio_edge_pi.x+t2*v[1], vio_edge_pi.y+t2*v[2])

	# in rotation left, only need to check with t1
	if t1 > 0 && t1 < 1
		# this is the one to use
		# rotation to right is good
		return angle_between_rotation(point, myPoint2d(vio_edge_pi.x, vio_edge_pi.y), new_point)
	elseif t1 > 1
		# need to check with previous edge
		return rotate_left_against_edge(peg, socket, current, point, vio[length(vio)].point, vio[length(vio)].edge+1)
	else
		# other cases,
		# cannot work?
		# should not happen
		return 0
	end

	return 0

end


function rotate_left(peg::myPeg, socket::mySocket, vio::Array{Int}, current::Int, point::myPoint2d)
	vio_point = peg.points[vio[length(vio)].point]
	vio_edge_id = socket.edges[vio[length(vio)].edge]

	vio_edge_pi = socket.points[vio_edge_id]
	vio_edge_pj = socket.points[vio_edge_id]

	# now that we have got all the points
	# compute how far the vio_point is from the rotation center
	dist = distance(point, vio_point)

	# rotate until contact with edge
	# then detect if on the edge, or violate another edge

	t1, t2 = line_circle_intersection(vio_edge_pi, vio_edge_pj, point, dist)

	# because rotating right, should only care the one
	# that is between 0 and 1
	# if both would be between 0 and 1
	# use the smaller one
	# from the function
	# t1 > t2 strictly!

	v = [vio_edge_pj.x - vio_edge_pi.x, vio_edge_pj.y - vio_edge_pi.y]
	new_point = myPoint2d(vio_edge_pi.x+t2*v[1], vio_edge_pi.y+t2*v[2])

	# in rotation left, only need to check with t1
	if t1 > 0 && t1 < 1
		# this is the one to use
		# rotation to right is good
		return angle_between_rotation(point, myPoint2d(vio_edge_pi.x, vio_edge_pi.y), new_point)
	elseif t1 > 1
		# need to check with previous edge
		return rotate_left_against_edge(peg, socket, current, point, vio[length(vio)].point, vio[length(vio)].edge+1)
	else
		# other cases,
		# cannot work?
		# should not happen
		return 0
	end

	return 0
end

function rotate_right_against_edge(peg::myPeg, socket::mySocket, current::Int, point::myPoint2d,
		contact_point::Int, contact_edge::Int)


	vio_point = peg.points[contact_point]
	vio_edge_id = socket.edges[contact_edge]

	vio_edge_pi = socket.points[vio_edge_id]
	vio_edge_pj = socket.points[vio_edge_id]

	dist = distance(point, vio_point)

	t1, t2 = line_circle_intersection(vio_edge_pi, vio_edge_pj, point, dist)

	v = [vio_edge_pj.x - vio_edge_pi.x, vio_edge_pj.y - vio_edge_pi.y]
	new_point = myPoint2d(vio_edge_pi.x+t2*v[1], vio_edge_pi.y+t2*v[2])

	# in rotation right, only need to check with t2
	if t2 > 0 && t2 < 1
		# this is the one to use
		# rotation to right is good
		return angle_between_rotation(point, myPoint2d(vio_edge_pi.x, vio_edge_pi.y), new_point)
	elseif t2 < 0
		# need to check with previous edge
		return rotate_right_against_edge(peg, socket, current, point, vio[length(vio)].point, vio[length(vio)].edge-1)
	else
		# other cases,
		# cannot work?
		# should not happen
		return 0
	end

	return 0

end


function rotate_right(peg::myPeg, socket::mySocket, vio::Array{Int}, current::Int, point::myPoint2d)
	# first, find the closest violation
	# then, make sure it is valid
	# then, check if other are in violation

	# because this is rotating right
	# the closest point in violation is
	# at the end of the violation list

	#

	# point: the rotation center

	# first, find the closest point
	# closest = vio[len(vio)].point # this is the id
	vio_point = peg.points[vio[length(vio)].point]
	vio_edge_id = socket.edges[vio[length(vio)].edge]

	vio_edge_pi = socket.points[vio_edge_id]
	vio_edge_pj = socket.points[vio_edge_id]

	# now that we have got all the points
	# compute how far the vio_point is from the rotation center
	dist = distance(point, vio_point)

	# rotate until contact with edge
	# then detect if on the edge, or violate another edge

	t1, t2 = line_circle_intersection(vio_edge_pi, vio_edge_pj, point, dist)

	# because rotating right, should only care the one
	# that is between 0 and 1
	# if both would be between 0 and 1
	# use the smaller one
	# from the function
	# t1 > t2 strictly!

	v = [vio_edge_pj.x - vio_edge_pi.x, vio_edge_pj.y - vio_edge_pi.y]
	new_point = myPoint2d(vio_edge_pi.x+t2*v[1], vio_edge_pi.y+t2*v[2])

	# in rotation right, only need to check with t2
	if t2 > 0 && t2 < 1
		# this is the one to use
		# rotation to right is good
		return angle_between_rotation(point, myPoint2d(vio_edge_pi.x, vio_edge_pi.y), new_point)
	elseif t2 < 0
		# need to check with previous edge
		return rotate_right_against_edge(peg, socket, current, point, vio[length(vio)].point, vio[length(vio)].edge-1)
	else
		# other cases,
		# cannot work?
		# should not happen
		return 0
	end

	return 0


end


function rotate_to_valid_config(peg::myPeg, socket::mySocket, vio::Array{Int}, current::Int,
		point::myPoint2d)
	# see how to rotate


	# next, to check if there's any conflict any more
	# after the rotation

	if vio[length(vio)].point < current
		# left
		result = rotate_left(peg, socket, vio, current, point)
	elseif vio[1].point > current
		# right
		result = rotate_right(peg, socket, vio, current, point)
		# because rotating to right, angle is negative
		result = - result

		# need to check if there's any other conflict
	end


	# maybe rotate the peg by that angle
	# then, get new configuration, test again
	config.th = result
	left, right, vio = check_violation_in_seq(peg, config, socket, current)
	if !left && !right
		# no collision
		return true, config
	elseif !left || !right
		# need further rotation
		return rotate_to_valid_config(peg, socket, vio, current, myPoint2d(v[1], v[2]))
	end

	return false, nothing

end


function find_valid_for_contact_pair(peg::myPeg, socket::mySocket, c_p::contact_pair)

	# the idea is to first test contact on both end of the edge;
	# and if either exist a valid rotation, then return true;
	# if both are invalid, search in a binary fashion, until find a region / location on edge
	# so that the peg can be collision free
	# or conflict: collision on both sides of the contact point
	# so rotation cannot work
	#

	# can be viewed as an alternating optimization
	# to find valid configuration for a peg

	current_cp = c_p.c_p.index
	current_ce = c_p.c_e.index

	# first, no rotation, only translation, on both end of the edge

	# contact point, in peg frame
	cx = c_p.c_p.x
	cy = c_p.c_p.y

	# contact edge starting point
	pix = c_p.c_e.s.x
	piy = c_p.c_e.s.y

	# contact edge ending point
	pjx = c_p.c_e.e.x
	pjy = c_p.c_e.e.y


	# first, translate to pi, and test violation

	v = [pix-cx, piy-cy] # the vector to translate contact point to pi
	config = configSE2(v[1], v[2], 0)
	# test the violation
	left, right, vios = check_violation_in_seq(peg, config, socket, current_cp)

	if !left && !right
		# both have no collisions
		return true
	end

	# if left && right
	# 	# both have collisions
	# 	# then, no rotation can fix this at this point
	# 	# move on to next text?
	# 	# so this should not be an if statement
	# end

	# if either only left or only right,
	# see if rotation can fix it
	if !left || !right
		# only one of them is false
		# and only one of them is false
		# see if rotation can find a possible configuration
		result, config = rotate_to_valid_config(peg, socket, vio, current, myPoint2d(v[1], v[2]))
	end


	if result
		# exist configuration
		return true, config
	end






	# then, translate to pj, and test violation

	v = [pjx-cx, pjy-cy] # the vector to translate contact point to pi
	config = configSE2(v[1], v[2], 0)
	# test the violation
	left, right, vios = check_violation_in_seq(peg, config, socket, current_cp)

	if !left && !right
		# both have no collisions
		return true
	end

	# if left && right
	# 	# both have collisions
	# 	# then, no rotation can fix this at this point
	# 	# move on to next text?
	# 	# so this should not be an if statement
	# end

	# if either only left or only right,
	# see if rotation can fix it
	if !left || !right
		# only one of them is false
		# and only one of them is false
		# see if rotation can find a possible configuration
		result, config = rotate_to_valid_config(peg, socket, vio, current, myPoint2d(v[1], v[2]))
	end


	if result
		# exist configuration
		return true, config
	end


	# above test will finish if there exist valid configuration
	# at either endpoint of the edge; if no, then we know this contact-mode
	# is not fully free along the entire edge
	# therefore, need to recursive find possible
	# intermediate point along the edge, that is collision free


	# last, translate to mid between pi and pj, test violation
	#


	# if valid, done, return;


	# if not valid, test which region is possible, and do a recursive search


end


function inward_normal(point1::myType.myPoint2d, point2::myType.myPoint2d)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function inward_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function compute_signed_angle(v1, v2)

	vn = [0,0,1]

	#os = LinearAlgebra.dot(v1,v2)/v1[1]^2+v2

	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))

	return th


end

function test_cp_ag()

	v1 = [1,0,0]
	v2 = [-1,1,0]
	th = compute_signed_angle(v1,v2)
	println(th/pi)
end


function inside_poly(point::Indexed_point, p_list::Array{Indexed_point})
	th = 0
	x = point.x
	y = point.y
	n = length(p_list)

	for i = 1:n
		next = (i+1)%n
		if next == 0
			next = n
		end
		v1 = [p_list[i].x-x, p_list[i].y-y, 0]
		v2 = [p_list[next].x-x, p_list[next].y-y, 0]

		th += compute_signed_angle(v1, v2)
	end

	n = th/(2*pi)
	println(th,",",n)
	if all_close(n%2, 1)
		return true
	else
		return false
	end
end

function intersection_between_lines(l1::Array, l2::Array)
	if size(l1, 1) != size(l2, 1)
		println("Error! wrong dimension of lines. ")
		exit(1)
	end

	A = Matrix{Float64}(undef, 2, 2)
	A[1, 1] = l1[1]
	A[1, 2] = l1[2]
	A[2, 1] = l2[1]
	A[2, 2] = l2[2]
	B = Array{Float64}(undef, 2)
	B[1] = -l1[3]
	B[2] = -l2[3]

	result = LinearAlgebra.pinv(A) * B

	return result

end


function line_seg_inter(p_list1::Array{Indexed_point},p_list2::Array{Indexed_point})
	l1 = get_line_eq(p_list1[1],p_list1[2])
	l2 = get_line_eq(p_list2[2],p_list2[2])

	min_x1 = min(p_list1[1].x, p_list1[2].x)
	max_x1 = max(p_list1[1].x, p_list1[2].x)
	min_x2 = min(p_list2[1].x, p_list2[2].x)
	max_x2 = max(p_list2[1].x, p_list2[2].x)

	min_x = max(min_x1, min_x2)
	max_x = min(max_x1, max_x2)

	min_y1 = min(p_list1[1].y, p_list1[2].y)
	max_y1 = max(p_list1[1].y, p_list1[2].y)
	min_y2 = min(p_list2[1].y, p_list2[2].y)
	max_y2 = max(p_list2[1].y, p_list2[2].y)

	min_y = max(min_y1, min_y2)
	max_y = min(max_y1, max_y2)


	p = intersection_between_lines(l1,l2)

	if p[1]>=min_x && p[1]<=max_x && p[2]>=min_y &&p[2]<=max_y
		return 1
	else
		return 0
	end
end

function check_inside_socket(point::Indexed_point, socket::mySocket)
	#points are anticlockwise
	ps = Array{Indexed_point}(undef, 0)
	for i = 1:size(socket.points,1)
		push!(ps, socket.points[i])
	end

	return inside_poly(point, ps)

end

function check_outside_peg(point::Indexed_point, peg::myPeg)

	ps = Array{Indexed_point}(undef, 0)
	for i = 1:size(peg.points,1)
		push!(ps, peg.points[i])
	end

	return !inside_poly(point, ps)

end



function check_socket_violation(peg::myPeg, config::configSE2, socket::mySocket, current::Int)

	left = false
	right = false

	new_peg = Kinematics.get_peg_from_config(peg, config)
	vio = Array{c_pair}(undef, 0)
	n = size(socket.points, 1)

	for i = 1:n
		result = check_outside_peg(socket.points[i], new_peg)

		if !result

			#=
			last = i-1
			if last == 0
				last = n
			end

			l2 = [socket.points[socket.edges[last].s], socket.points[socket.edges[last].e]]
			=#
			l1 = [socket.points[socket.edges[i].s], socket.points[socket.edges[i].e]]


			for j = 1:size(peg.edges, 1)
				l2 = [peg.points[peg.edges[j].s], peg.points[peg.edges[j].e]]
				#=
				count = line_seg_inter(l1, l3) + line_seg_inter(l2, l3)
				=#
				if line_seg_inter(l1, l2) == 1
					push!(vio, c_pair(i, j, 0))

					if j < current
						left = true
					else
						right = true
					end

					break
				end
			end
		end

	end


	return left, right, vio

end

function check_peg_violation(peg::myPeg, config::configSE2, socket::mySocket, current::Int)
	left = false
	right = false

	ps = Kinematics.compute_points_from_config(peg, config)
	vio = Array{c_pair}(undef, 0)
	n = size(ps, 1)

	for i = 1:n
		if i == current
			continue
		end

		result = check_inside_socket(ps[i], socket)
		println(i,",",result)
		if !result

			l1 = [peg.points[peg.edges[i].s], peg.points[peg.edges[i].e]]


			for j = 1:size(socket.edges, 1)
				l2 = [peg.points[peg.edges[j].s], peg.points[peg.edges[j].e]]

				if line_seg_inter(l1, l2) == 1
					push!(vio, c_pair(i, j, 1))

					if i < current
						left = true
					else
						right = true
					end

					break
				end
			end
		end

	end


	return left, right, vio

end

function check_violation(peg::myPeg, config::configSE2, socket::mySocket, current::Int)
	left1, right1, vio1 = check_peg_violation(peg, config, socket, current)
	left2, right2, vio2 = check_socket_violation(peg, config, socket, current)


	vio = vcat(vio1,vio2)

	left = left1||left2
	right = right1||right2
	return left, right, vio
end

function mytest()
	peg,socket = Construct_Peg_Socket.read_from_file(1)
	points1=Array{myType.myPoint2d}(undef, 0)
	points2=Array{myType.myPoint2d}(undef, 0)

	for i = 1:size(peg.points,1)
		push!(points1, myType.myPoint2d(peg.points[i].x, peg.points[i].y))
	end
	for i = 1:size(socket.points,1)
		push!(points2, myType.myPoint2d(socket.points[i].x, socket.points[i].y))
	end



	Visualization.draw_polygon(points1)
	Visualization.draw_polygon(points2)

	pl1 = Array{myType.Indexed_point}(undef,2)
	pl2 = Array{myType.Indexed_point}(undef,2)

	pl1[1]=myType.Indexed_point(1, 0., 0.)
	pl1[2]=myType.Indexed_point(2, 1., 0.)
	pl2[1]=myType.Indexed_point(1, 2., 3.)
	pl2[2]=myType.Indexed_point(2, 1., 1.)


	#result = line_seg_inter(pl1,pl2)
	#config = configSE2(0.,0.,0.)
	#left, right, vio = check_violation(peg, config, socket, 1)

	#println(left, right, vio)

end

function test_cp(i::Int)
	peg, socket = Construct_Peg_Socket.read_from_file(i)

	contact_point = myType.Indexed_point(peg.points[3])
	contact_edge = myType.Indexed_edge(2, socket.points[2], socket.points[3])

	cp1 = myType.contact_pair(1, contact_point, contact_edge, 1)

	vx, vy, vth = find_valid_cp(peg, socket, cp1)


	# if I test with vx, vy, vth

	# configuration

	ps = Kinematics.compute_points_from_config(peg, myType.configSE2(vx, vy, vth))

	# ps = Kinematics.compute_points_from_config(peg, myType.configSE2(0, 0, 0))

	Visualization.display_peg_socket(socket, ps)

end

function signed_angle(th,x,y,p_x,p_y,p_x_next,p_y_next,pi_x,pi_y)
	#ps,pi_x,pi_y are parameters


		v1=[(cos(th) * p_x - sin(th) * p_y + x)-pi_x,
		(sin(th) *p_x + cos(th) * p_y + y)-pi_y]

		v2=[(cos(th) * p_x_next - sin(th) * p_y_next + x)-pi_x,
		(sin(th) *p_x_next + cos(th) * p_y_next + y)-pi_y]

		vn = [0,0,1]


		#os = LinearAlgebra.dot(v1,v2)/v1[1]^2+v2

		#tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


		th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))


		return th

end

function test_NLPV(x,th,p)

	return tan(p*(x - th))



end

function test_NL()

	m = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))
	JuMP.@variable(m, x)
	JuMP.@variable(m, th)

	a = [1,2,3,4]
	n=4
	JuMP.@NLparameter(m,p[i=1:4]==a[i])
	JuMP.register(m, :test, 3, test_NLPV, autodiff=true)

	JuMP.@NLexpression(m, exp1, sum(atan(test(x,th,p[i])) for i = 1:n))

	JuMP.@NLconstraint(m, abs(exp1-pi)<=0.0001)
	JuMP.@NLobjective(m, Min, x^2+th^2)
	status = JuMP.solve(m)

	vx = JuMP.getvalue(x)
	vth = JuMP.getvalue(th)

	println(vx,vth)



end



function test_cons(i)

	peg, socket = Construct_Peg_Socket.read_from_file(i)

	contact_point = myType.Indexed_point(peg.points[3])
	contact_edge = myType.Indexed_edge(2, socket.points[2], socket.points[3])

	c_p = myType.contact_pair(1, contact_point, contact_edge, 1)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))
	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)

	ps_x = socket.points[1].x
	ps_y = socket.points[1].y


	ps = Array{Indexed_point}(undef, 0)
	ps_next = Array{Indexed_point}(undef, 0)
	n = length(peg.points)
	for i = 1:n
		next = (1+i)%n
		if next==0
			next = n
		end
		push!(ps, peg.points[i])
		push!(ps_next, peg.points[next])

	end



	JuMP.@NLparameter(model, p_x[i=1:n]==ps[i].x)
	JuMP.@NLparameter(model, p_y[i=1:n]==ps[i].y)
	JuMP.@NLparameter(model, p_x_next[i=1:n]==ps_next[i].x)
	JuMP.@NLparameter(model, p_y_next[i=1:n]==ps_next[i].y)
	JuMP.@NLparameter(model, pi_x==ps_x)
	JuMP.@NLparameter(model, pi_y==ps_y)


	JuMP.register(model, :signed_angle, 9, signed_angle, autodiff=true)

	JuMP.@NLexpression(model, exp, sum(signed_angle(th,x,y,p_x[i],p_y[i],p_x_next[i],p_y_next[i],pi_x,pi_y)
	for i =1:n ))


	if c_p.mode==0 && i == c_p.c_p.index

		JuMP.@NLconstraint(model, abs(exp-pi)<=0.00001)
	else

		JuMP.@NLconstraint(model, abs(exp)<=0.00001)
	end


	JuMP.@NLobjective(model, Min, x*x+y*y+th*th)


	status = JuMP.solve(model)

	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	if status == :Optimal
		println(vx, vy, vth)
	else
		println(error)
	end




end


function find_valid_cp(peg::myPeg, socket::mySocket, c_p::contact_pair)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))
	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)

	for i = 1:length(socket.points)

		pi_x = socket.points[i].x
		pi_y = socket.points[i].y
		total = 0

		p_s = Array{Indexed_point}(undef, 0)


		for j = length(peg.points)

			#ps = peg.points[peg.edges[j].s]
			#pe = peg.points[peg.edges[j].e]

			#v_s = [(cos(th) * ps.x - sin(th) * ps.y + x)-pi_x,(sin(th) *ps.x + cos(th) * ps.y + y)-pi_y]
			#v_e = [(cos(th) * pe.x - sin(th) * pe.y + x)-pi_x,(sin(th) *pe.x + cos(th) * pe.y + y)-pi_y]

			#total += compute_signed_angle(v_s, v_e)
			push!(p_s, peg.points[j])


		end

		JuMP.@NLparameter(model, ps[i=1:length(p_s)]==p_s[i])

		JuMP.register(model, :sum_signed_angle, 6, sum_signed_angle, autodiff=true)

		JuMP.@NLconstraint(model, sum_signed_angle(th,x,y,ps,pi_x,pi_y)<=0.00001)


		if c_p.mode==0 && i == c_p.c_p.index

			JuMP.@NLconstraint(model, abs(sum_signed_angle(th,x,y,ps,pi_x,pi_y)-pi)<=0.00001)
		else

			JuMP.@NLconstraint(model, abs(sum_signed_angle(th,x,y,ps,pi_x,pi_y))<=0.00001)
		end
	end




	for i = 1:length(peg.points)

		pi_x = peg.points[i].x
		pi_y = peg.points[i].y
		total = 0
		n = length(socket.points)
		for j = n
			next = (j+1)%n
			if next ==0
				next = n
			end

			ps = socket.points[j]
			pe = socket.points[next]

			v_s = [(cos(th) * ps.x - sin(th) * ps.y + x)-pi_x,(sin(th) *ps.x + cos(th) * ps.y + y)-pi_y]
			v_e = [(cos(th) * pe.x - sin(th) * pe.y + x)-pi_x,(sin(th) *pe.x + cos(th) * pe.y + y)-pi_y]

			total += compute_signed_angle(v_s, v_e)
		end

		if c_p.model==1 && i == c_p.c_p.index

			JuMP.@NLconstraint(model,
				abs(total - pi) <= 0.00001)
		else
			JuMP.@NLconstraint(model,
				abs(total - 2*pi) <= 0.00001)

		end
	end

	JuMP.@NLobjective(model, Min, th*th)


	status = JuMP.solve(model)

	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	if status == :Optimal
		return vx, vy, vth
	end

	return nothing, nothing, nothing


end

function find_valid_pair_constraints(peg::myPeg, socket::mySocket, c_p::contact_pair)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	# unknown variables: the configuration
	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)

	for i = 1:length(socket.edges)
		p_i = socket.points[socket.edges[i].s]
		# get the inward normal
		n = inward_normal(p_i, socket.points[socket.edges[i].e])

		for j = 1:length(peg.points)
			if j == c_p.c_p.index
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			JuMP.@NLconstraint(model,
				n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
				n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0)

		end
	end



	p_i = c_p.c_e.s
	n = inward_normal(p_i, c_p.c_e.e)

	p = c_p.c_p
	# last add the equality constraint of the contact pair
	JuMP.@NLconstraint(model,
		n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
		n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) == 0);





	# minimize the absolute value of theta

	JuMP.@NLobjective(model, Min, th*th)
	# JuMP.@variable(model, z)
	# JuMP.@constraint(model, z >= th)
	# JuMP.@constraint(model, z <= -th)

	# JuMP.@objective(model, Min, z)

	status = JuMP.solve(model)

	# println
	# println("got ", JuMP.getobjectivevalue(model), " at ",
	# 	JuMP.getvalue(x), ", ", JuMP.getvalue(y), ", ", JuMP.getvalue(th))
	# end print
	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	vio = false

	# println("status: ", status)
	# println("values: ", vx, ", ", vy, ", ", vth)
	# println("test with all 0")

	# vx = 0
	# vy = 0
	# vth = 0


	for i = 1:length(socket.edges)
		p_i = socket.points[socket.edges[i].s]
		# get the inward normal
		n = inward_normal(p_i, socket.points[socket.edges[i].e])

		for j = 1:length(peg.points)
			if j == c_p.c_p.index
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			# println(n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
			# 	n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) )
			if n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
				n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) < -0.01
				vio = true
			end

		end
	end

	p_i = c_p.c_e.s
	n = inward_normal(p_i, c_p.c_e.e)

	p = c_p.c_p
	# last add the equality constraint of the contact pair
	# println("eq: ", n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
	# 	n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y) );
	if !all_close(n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
			n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y), 0)
		vio = true
	end




	if status == :Optimal && !vio
		return vx, vy, vth
	end

	return nothing, nothing, nothing

end

function point_in_contact_mode(index::Int, c_m::contact_mode)
	for i = 1:length(c_m.cps)

		if index == c_m.cps[i].c_p.index
			return true
		end
	end

	return false

end

function all_close(a, b)
	if abs(a - b) < 0.001
		return true
	end
	return false;

end



function find_valid_contact_mode(peg::myPeg, socket::mySocket, c_m::contact_mode)
	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.001, maxtime=5))

	# unknown variables: the configuration
	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)



	for i = 1:length(socket.edges)
		p_i = socket.points[socket.edges[i].s]
		# get the inward normal
		n = inward_normal(p_i, socket.points[socket.edges[i].e])

		for j = 1:length(peg.points)
			if point_in_contact_mode(j, c_m)
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			JuMP.@NLconstraint(model,
				n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
				n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0)

		end
	end


	# NEED TO MODIFy
	for i = 1:length(c_m.cps)
		p_i = c_m.cps[i].c_e.s
		n = inward_normal(p_i, c_m.cps[i].c_e.e)

		p = c_m.cps[i].c_p
		# last add the equality constraint of the contact pair
		# JuMP.@NLconstraint(model,
		# 	n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
		# 	n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) == 0);
		JuMP.@NLconstraint(model,
			n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
			n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) >= 0);
		JuMP.@NLconstraint(model,
			n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
			n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) <= 0);
	end

	# minimize the absolute value of theta

	JuMP.@NLobjective(model, Min, th*th + x*x + y*y)
	# JuMP.@objective(model, Min, 100)

	status = JuMP.solve(model)

	# println
	# println("got ", JuMP.getobjectivevalue(model), " at ",
	# 	JuMP.getvalue(x), ", ", JuMP.getvalue(y), ", ", JuMP.getvalue(th))
	# end print

	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	vio = false

	# ########### test constraints
	for i = 1:length(socket.edges)
		p_i = socket.points[socket.edges[i].s]
		# get the inward normal
		n = inward_normal(p_i, socket.points[socket.edges[i].e])

		for j = 1:length(peg.points)

			if point_in_contact_mode(j, c_m)
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			if n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
				n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) < -0.01
				# println("got here non negative? ")
				# println(n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
				# n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y))
				vio = true
			end

		end
	end


	for i = 1:length(c_m.cps)
		p_i = c_m.cps[i].c_e.s
		n = inward_normal(p_i, c_m.cps[i].c_e.e)

		p = c_m.cps[i].c_p
		# last add the equality constraint of the contact pair
		if !all_close(n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
				n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y), 0)
			# println("got here inq")
			# println(n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
			# 	n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y))
			vio = true

		end
	end



	################## end test of constraints




	if status == :Optimal && !vio
		return vx, vy, vth
	end

	return nothing, nothing, nothing

	# if status == :Optimal
	# 	return true
	# end

	# return false;
end

function test_contact_pair(i::Int)
	peg, socket = Construct_Peg_Socket.read_from_file(i)

	contact_point = myType.Indexed_point(peg.points[3])
	contact_edge = myType.Indexed_edge(2, socket.points[2], socket.points[3])

	cp1 = myType.contact_pair(1, contact_point, contact_edge)

	vx, vy, vth = find_valid_pair_constraints(peg, socket, cp1)


	# if I test with vx, vy, vth

	# configuration

	ps = Kinematics.compute_points_from_config(peg, myType.configSE2(vx, vy, vth))

	# ps = Kinematics.compute_points_from_config(peg, myType.configSE2(0, 0, 0))

	Visualization.display_peg_socket(socket, ps)

end

function test_constraint(c_i::Int, e_j::Int)
	# this is a hard coded peg, for testing only
	peg, socket = Construct_Peg_Socket.read_from_file(1)

	println("contact point: ", peg.points[c_i].x, ", ", peg.points[c_i].y)
	println("contact edge, from ", socket.points[socket.edges[e_j].s], " to ",
		socket.points[socket.edges[e_j].e])

	p_i = socket.points[socket.edges[e_j].s]
	p_j = socket.points[socket.edges[e_j].e]

	p = peg.points[c_i]

	# now we have got the points

	n = inward_normal(p_i, p_j)

	# now test the relation

	x = 0
	y = 0
	th = 0

	println(n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
				n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) )


end

function test_contact_mode(cis::Array{Int}, ejs::Array{Int})
	# this is a hard coded peg, for testing only
	peg, socket = Construct_Peg_Socket.read_from_file(1)

	ps = []
	es = []

	# contact_point = myType.Indexed_point(peg.points[3])
	# contact_edge = myType.Indexed_edge(2, socket.points[2], socket.points[3])

	# insert points

	for i = 1:size(cis, 1)
		push!(ps, peg.points[cis[i]])
		push!(es, myType.Indexed_edge(ejs[i],
			socket.points[socket.edges[ejs[i]].s], socket.points[socket.edges[ejs[i]].e]))
	end

	# now construct pairs

	cps = Array{myType.contact_pair}(undef, 0)
	for i = 1:size(cis, 1)
		push!(cps, myType.contact_pair(i, ps[i], es[i]))
	end

	c_m = myType.contact_mode(cps)





	vx, vy, vth = find_valid_contact_mode(peg, socket, c_m)

	if vx == nothing
		println("Error! not valid")
		return
	end

	ps = Kinematics.compute_points_from_config(peg, myType.configSE2(vx, vy, vth))

	# ps = Kinematics.compute_points_from_config(peg, myType.configSE2(0, 0, 0))

	Visualization.display_peg_socket(socket, ps)


end




function find_all_pairs(peg::myType.myPeg, socket::myType.mySocket)
	# for this, should select and test all the possible pairs

	lp = length(peg.points)
	le = length(socket.edges)
	se = myType.mySocketEdges(socket)

	ps = []

	for i = 1:lp
		for j = 1:le
			# try each possible combination
			# test pair

			c_p = myType.contact_pair((i-1)*lp + j, peg.points[i], se.edges[j])
			vx, vy, vth = find_valid_pair_constraints(peg, socket, c_p)
			if vx == nothing
				continue
			end


			# pts = Kinematics.compute_points_from_config(peg, myType.configSE2(vx, vy, vth))


			# Visualization.display_peg_socket(socket, pts)
			# sleep(1)

			# otherwise, it exists,
			push!(ps, [i, j])
		end
	end

	return ps

	# println(ps)

end

# assumption: peg is convex
# i.e. all the contact points forms a convex hull of the peg
# and all contact points are on the convex hull

function validate_new_pair(cps, current, newcp)

	# current can be a list of pairs
	# new is the new pair

	if length(current) == 1
		# do not select the second pair that
		# share the same contact point, but not adjacent edges
		cpi = cps[current[1]][1]
		cei = cps[current[1]][2]

		cpj = cps[newcp][1]
		cej = cps[newcp][2]

		if cpi == cpj && abs(cei - cej) != 1
			# same contact point, not adjacent edge
			return false
		end

		if cei == cej && abs(cpi - cpj) != 1
			# same contact edge
			# so, non-adjacent contact points cannot contact same edge
			return false
		end

		# any other possibility for sanity check?
		return true
	end

	# if current is longer than 1
	# need to check that the new pair is not colliding with any of the current pairs

	c_p = cps[newcp][1]
	c_e = cps[newcp][2]

	appearance = Dict{Int, Int}()

	for i = 1:length(current)
		cpi = cps[current[i]][1]
		cei = cps[current[i]][2]

		if haskey(appearance, cpi)
			appearance[cpi] += 1
		else
			appearance[cpi] = 1
		end

		if cpi == c_p && abs(cei - c_e) != 1
			# same contact point, not adjacent edge
			return false
		end

		if cpi == c_p && appearance[cpi] > 1
			return false
		end

		if cei == c_e && abs(cpi - c_p) != 1
			# same contact edge
			# so, non-adjacent contact points cannot contact same edge
			return false
		end




	end






	return true

end


function pick_k_from_array(k::Int, indicies::Array{Int}, result::Array{Int}, total, cps)
	if k == 0
		push!(total, result)
	end

	# total = []

	for i = 1:length(indicies)
		if indicies[i] < result[length(result)]
			continue
		end


		if ! validate_new_pair(cps, result, indicies[i])
			continue
		end

		nids = deepcopy(indicies)
		nrs = deepcopy(result)

		# println("before add: ", nrs)
		push!(nrs, nids[i])
		# println("after add: ", nrs)


		filter!(x->x≠nids[i], nids)
		# println("nids after filter is: ", nids)

		outcome = pick_k_from_array(k-1, nids, nrs, total, cps)
		if outcome == nothing
			continue
		elseif length(outcome) == 0
			continue
		end
		# println("temp: ", outcome)
		# push!(total, outcome)
	end

	return total

end


function pick_k(k::Int, cps)
	# l is the list
	# pick k numbers from l

	indicies = Array{Int}(undef, 0)
	for i = 1:length(cps)
		push!(indicies, i)
	end

	# now, pick k of them
	# result = Array{Int}(undef, 0)

	total = []

	for i = 1:length(indicies)
		result = Array{Int}(undef, 0)
		push!(result, i)
		nids = deepcopy(indicies)

		filter!(x->x≠i, nids)

		pick_k_from_array(k-1, nids, result, total, cps)

	end


	return total

end




function validate_mode(list, peg::myType.myPeg, socket::myType.mySocket, c_p_s)

	se = myType.mySocketEdges(socket)
	# println("enter: ", length(list))
	# list is a mode
	# not a list of modes

	invalid = []

	for mode in list
		# mode, map the mode

		contact_pairs = Array{myType.contact_pair}(undef, 0)
		count = 1
		for i in mode


			push!(contact_pairs, myType.contact_pair(count, peg.points[c_p_s[i][1]], se.edges[c_p_s[i][2]]))

			count += 1

		end

		c_m = myType.contact_mode(contact_pairs)
		# println("mode: ", mode)
		# now we have the mode
		vx, vy, vth = find_valid_contact_mode(peg, socket, c_m)
		if vx == nothing
			# remove this from the list
			push!(invalid, mode)
		end

	end

	for mode in invalid
		filter!(x->x≠mode, list)
	end

	# println("exiting: ", length(list))

end



function all_combinations(peg::myType.myPeg, socket::myType.mySocket)
	# for each contact point, form with each of the edge
	# also, all possible number of pairs
	# how to loop?

	# first loop over the possible number of pairs

	# contact pairs
	c_p_s = find_all_pairs(peg, socket)
	# println(length(c_p_s))
	println(c_p_s)
	# println("length: ", length(peg.points))
	# for all possible combinations


	for i = 2:length(c_p_s)
		# the most possible is that all contact points are
		# in contact with one of the edge

		# pick i pairs, from the result of previous set and
		# problem: a single contact point can contact 2 edges at the same time
		# but, in fact that should be a special situation, as that
		# either is a transition point, or a sink
		# println("i: ", i)
		total = pick_k(i, c_p_s)
		# println(length(total))
		# println(total)
		# break

		validate_mode(total, peg, socket, c_p_s)

	end

	# for every result in the mode
	# test if the mode is possible

	# total = pick_k(5, c_p_s)

end


function adjacent_modes(m1::myType.contact_mode, m2::myType.contact_mode)

	l1 = length(m1.cps)
	l2 = length(m2.cps)

	# if too many change of pairs, not adjacent
	# it is possible that when a contact point
	# switches contact edge
	# two new contact points contact new edges
	# so, the following criteria is not totally correct


	# if abs(l1 - l2) > 1
	# 	return false
	# end

	# if same length
	#



end


function sort_contact_mode(cm::myType.contact_mode)

	# the goal is to sort the contact_mode
	# based on the increasing order of contact points
	# and contact edges;
	# first, contact points are sorted based on increasing orders;
	# and if the contact points are the same
	# the contact edges are sorted based on increasing order.

	# can use insertion sort
	# or merge sort
	# merge sort can be fast
	# but most of the contact mode is short
	# so, does not make much difference


	l = length(cm.cps)

	# first approach: insertion sort
	for i = 2:l
		# The pair to be moved
		temp = cm.cps[i]

		for j = i-1:-1:1
			#
			if (cm.cps[j].c_p.index > temp.c_p.index)
				cm.cps[j+1] = cm.cps[j]
			else
				# loop over contact points with same id
				# sort contact edges

				for k = j:-1:1
					if (cm.cps[k].c_p.index == temp.c_p.index &&
							cm.cps[k].c_e.index > temp.c_e.index)
						cm.cps[k+1] = cm.cps[k]
					else
						cm.cps[k+1] = temp
						break
					end
				end
				break
			end
		end

	end

	return

end


function is_neighbor_modes(cm1::myType.contact_mode, cm2::myType.contact_mode)

	# off by length at most 1
	# cannot be all different
	# differ by at most one pair

	l1 = length(cm1.cps)
	l2 = length(cm2.cps)

	if abs(l1 - l2) > 1
		return false;
	end

	# else if the number of difference is larger than 1
	# the two modes cannot be neighbors

	sort_contact_mode(cm1)
	sort_contact_mode(cm2)

	# after sort
	# compare

	num_diff = 0
	for i = 1:l1
		current = cm1.cps[i]
		for j = 1:l2
			if cm2.cps[j].c_p.index < current.c_p.index
				continue
			elseif cm2.cps[j].c_p.index == current.c_p.index &&
				cm2.cps[j].c_e.index < current.c_e.index
				continue
			elseif cm2.cps[j].c_p.index == current.c_p.index &&
				cm2.cps[j].c_e.index == current.c_e.index
				break
			elseif cm2.cps[j].c_p.index == current.c_p.index &&
				cm2.cps[j].c_e.index > current.c_e.index
				num_diff += 1
				break
			elseif cm2.cps[j].c_p.index > current.c_p.index
				num_diff += 1
				break
			end
		end
	end

	if num_diff > 1
		return false
	end

	return true

end


# now that if we can say two modes are neighbors;
# if both modes are valid
# they can be connected
# now, compute distances

function dist_between_modes(cm1::myType.contact_mode, cm2::myType.contact_mode)
	# need to use force direction
	# the configuration
	# the peg
	# to determine the distance between to modes

end

function dist_cp_ce(c_p::myType.Indexed_point, c_e::myType.Indexed_edge,
					peg::myType.myPeg, socket::myType.mySocket)
	# cpw = [c_p.x * cos(config.th) - c_p.y * sin(config.th) + config.x,
	# 		c_p.x * sin(config.th) + c_p.y * cos(config.th) + config.y]

	cpw = [c_p.x, c_p.y]

	in_normal = inward_normal(c_e.s, c_e.e)

	v = [cpw[1] - c_e.s.x, cpw[2] - c_e.s.y]
	dist = v[1] * in_normal[1] + v[2] * in_normal[2]

	return abs(dist)
end


function rescale_vector(v::Array, target_scale=1)

	scale = 0
	for i = 1:length(v)
		scale += v[i] * v[i]
	end

	scale = sqrt(scale)

	factor = target_scale / scale

	u = deepcopy(v)
	for i = 1:length(v)
		u[i] = u[i] * factor
	end

	return u



end



# valid direction from cm1 to cm2
# based on the socket, peg, configuration
# and force direction
# now, for special consideration of the force
# if F = [0, 0, 0]: force can be any direction
# so, no need to test
# define special force direction
# being along the peg, with $\alpha$ error
# F = [0, 0, theta]
function valid_direction(cm1::myType.contact_mode, cm2::myType.contact_mode,
						peg::myType.myPeg, socket::myType.mySocket,
						config::myType.configSE2, F::Array)

	# if (F[1] == 0 && F[2] == 0 && F[3] == config.th)
	# 	# special force direction
	# 	# along peg direction
	# 	# with error angle
	# 	#
	# end

	# need to first consider the summation of forces
	# for cm1
	# then, compute the difference between cm1 and cm2
	# see if the difference can be caused by the summation of the force

	# first, find the different or additional pairs
	cm1_copy = deepcopy(cm1)
	cm2_copy = deepcopy(cm2)

	# now, loop over cm1, find the same pair in cm2
	# if found pair, remove both
	# if no match found after the entire search
	# the remaining pairs are the difference

	# should not use for loop
	# due to the fact that we will remove things

	i = 1
	while i <= length(cm1_copy.cps)
		j = 1
		current = cm1_copy.cps[i]
		found_match = false
		while j <= length(cm2_copy.cps)
			if current.c_p.index == cm2_copy.cps[j].c_p.index &&
				current.c_e.index == cm2_copy.cps[j].c_e.index
				# found matching pairs
				# remove both
				#
				# remove element in cm1 at i
				# remove element in cm2 at j
				deleteat!(cm1_copy.cps, i)
				deleteat!(cm2_copy.cps, j)
				found_match = true
				break
			end
			j += 1
		end
		if !found_match
			i += 1
		end
	end

	#

	# now, found the different pair
	# see whether towards that direction
	# the basic relation is a contact point towards an edge
	# so, based on the difference, first find the contact point and
	# contact edge
	contact_point = nothing
	contact_edge = nothing

	if length(cm1_copy.cps) == 0
		# cm1 has the extra
		contact_point = cm2_copy.cps[1].c_p
		contact_edge = cm2_copy.cps[1].c_e
	elseif length(cm2_copy.cps) == 0
		# cm2 has the extra
		contact_point = cm1_copy.cps[1].c_p
		contact_edge = cm1_copy.cps[1].c_e
	else
		# both have one,
		# find the difference
		# if the same is the contact point
		if cm1_copy.cps[1].c_p.index == cm2_copy.cps[1].c_p.index
			contact_point = cm1_copy.cps[1].c_p
			contact_edge = cm2_copy.cps[1].c_e
		else
			contact_point = cm2_copy.cps[1].c_p
			contact_edge = cm2_copy.cps[1].c_e
		end
	end

	# now that we have found the target contact point and contact edge
	# need to see if at every contact point
	# the force direction will reduce the distance between
	# the contact point and the contact edge
	# won't occur the case where
	# different contact point have different force direction
	# that will create torque
	# that cannot happen

	# first, compute the distance between a contact point
	# and a contact edge

	px = contact_point.x * cos(config.th) - contact_point.y * sin(config.th) + config.x
	py = contact_point.x * sin(config.th) + contact_point.y * cos(config.th) + config.y

	c_p = myType.Indexed_point(contact_point.index, px, py)

	dist = dist_cp_ce(c_p, contact_edge, peg, socket)
	for i = 1:length(cm1.cps)

		F_p = force_between_pair(cm1.cps[i], peg, socket, config, F)

		# this is the force at contact point
		#

		# need to normalize F_p so that it does not overshoot

		# this F_s is the scaled force
		F_s = rescale_vector(F, dist / 2)

		ncp = myType.Indexed_point(c_p.index, c_p.x + F_s[1],
								c_p.y + F_s[2])

		n_dist = dist_cp_ce(ncp, contact_edge, config, peg, socket)

		if (n_dist > dist)
			# increased distance
			return false
		end

	end


	return true


end

function force_between_pair(c_p::myType.contact_pair, peg::myType.myPeg,
				socket::myType.mySocket, config::myType.configSE2, F::Array)

	# compute the force between any contact pairs
	# given the orientation, force direction;
	# and the contact point location, contact edge angle
	# friction coefficient
	# one can compute the combined force direction
	# between a contact pair

	# force direction is along the given direction
	# as the block is pushed by a constant force

	# first need to compute the angle of the edge

	epi = c_p.c_e.s;
	epj = c_p.c_e.e;

	# direction vector
	edge_dir = [epj.x-epi.x, epj.y-epi.y]
	# inward normal
	i_n = inward_normal(epi, epj)
	o_n = negate_vector(i_n)

	# Now, compute the force
	# if it is a straight vector
	# or it is along the direction of the peg
	# in a cone

	if (F[1] == 0 && F[2] == 0 && F[3] == config.th)
		# special force direction
		# along peg direction
		# with error angle
		#

		# this is within a cone
		# need to make sure that the worst case is still ok

		# first compute the threshold force direction
		Ft = [(CONST.friction_coeff * o_n[2] - edge_dir[2]) / (edge_dir[1] - CONST.friction_coeff * o_n[1]), 1]

		ang = atan(Ft[2], Ft[1])

		ang_diff = angle_diff(ang, F[3])

		if ang_diff > CONST.force_cone
			# not in range
			# now, either in the direction of the edge
			# or in reverse

			ang = atan(edge_dir[2], edge_dir[1])

			a_d = angle_diff(ang, F[3]);
			if a_d > pi/2
				return negate_vector(edge_dir)
			else
				return edge_dir
			end

		else
			# in range
			# so can get stuck
			return [0, 0]
		end





	else
		# this is normal situation
		# force is a vector
		# can directly compute the projection
		# to different components
		# use the coefficient

		Fn = F[1] * o_n[1] + F[2] * o_n[2]

		Fe = F[1] * edge_dir[1] + F[2] * edge_dir[2]

		if Fe > CONST.friction_coeff * Fn
			# can move
			if Fe > 0
				return edge_dir
			else
				return negate_vector(edge_dir)
			end
		else
			return [0, 0]
		end

	end




end



function test(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)

	# find_all_pairs(peg, socket)
	all_combinations(peg, socket)

end


function test_sort()
	cps = []
	ces = []

	cp_s = []
	for i = 1:5
		push!(cps, myType.Indexed_point(i, 0.0, 0.0))
	end

	for i = 1:4
		push!(ces, myType.Indexed_edge(i, cps[i], cps[i+1]))
	end

	for i = 1:5
       for j = 1:4
           push!(cp_s, myType.contact_pair((i-1)*5+j-1, cps[i], ces[j]))
       end
   end

   a = [cp_s[1], cp_s[5], cp_s[2], cp_s[20], cp_s[18], cp_s[3], cp_s[15]]

   cm = myType.contact_mode(a)

   for i = 1:length(cm.cps)
	   println(cm.cps[i].c_p.index, ", ", cm.cps[i].c_e.index)
   end

   sort_contact_mode(cm)
   println("after sort")

   for i = 1:length(cm.cps)
	   println(cm.cps[i].c_p.index, ", ", cm.cps[i].c_e.index)
   end

end

function test_find_diff()

	# first generate some contact pairs and contact modes
	# make sure that the cm1 and cm2 is off at most 1

	cps = []
	ces = []

	cp_s = []
	for i = 1:5
		push!(cps, myType.Indexed_point(i, 0.0, 0.0))
	end

	for i = 1:4
		push!(ces, myType.Indexed_edge(i, cps[i], cps[i+1]))
	end

	for i = 1:5
       for j = 1:4
           push!(cp_s, myType.contact_pair((i-1)*5+j-1, cps[i], ces[j]))
       end
    end

    a = [cp_s[1], cp_s[2], cp_s[3], cp_s[4], cp_s[6]]

    b = [cp_s[1], cp_s[3], cp_s[4], cp_s[5], cp_s[6]]

    cm1 = myType.contact_mode(a)
    cm2 = myType.contact_mode(b)

    println(cp_s[2].id, ", ", cp_s[5].id)


	#
	#
	##################

	cm1_copy = deepcopy(cm1)
	cm2_copy = deepcopy(cm2)

	# now, loop over cm1, find the same pair in cm2
	# if found pair, remove both
	# if no match found after the entire search
	# the remaining pairs are the difference

	# should not use for loop
	# due to the fact that we will remove things

	i = 1
	while i <= length(cm1_copy.cps)
		j = 1
		current = cm1_copy.cps[i]
		found_match = false
		while j <= length(cm2_copy.cps)
			if current.c_p.index == cm2_copy.cps[j].c_p.index &&
				current.c_e.index == cm2_copy.cps[j].c_e.index
				# found matching pairs
				# remove both
				#
				# remove element in cm1 at i
				# remove element in cm2 at j
				deleteat!(cm1_copy.cps, i)
				deleteat!(cm2_copy.cps, j)
				found_match = true
				break
			end
			j += 1
		end
		if !found_match
			i += 1
		end
	end

	#
	println(cm1_copy.cps)
	#
	println(cm2_copy.cps)

	# loop over cm1

end


end #module
