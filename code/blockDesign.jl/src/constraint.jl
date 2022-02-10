module Constraints

import LinearAlgebra
import Luxor

import ..Construct_Peg_Socket
import ..Kinematics

using ..myType


# this module is used to check different constraints and their satisfaction status


# first batch of tests:
# whether all points on the peg is on one side of the socket
#
function all_close(a, b)
	if abs(a - b) < 0.0001
		return true
	end
	return false;

end


function get_line_eq(p1::Array, p2::Array)
    A = (p1[2] - p2[2])
    B = (p2[1] - p1[1])
    C = (p1[1]*p2[2] - p2[1]*p1[2])

    #line equation is Ax+By+c = 0
    return [A, B, C]
end

function get_line_eq(p1::Indexed_point, p2::Indexed_point)
    A = (p1.y - p2.y)
    B = (p2.x - p1.x)
    C = (p1.x*p2.y - p2.x*p1.y)
    return [A, B, C]
end

function get_two_line_intersection(L1::Array, L2::Array)
    D  = L1[1] * L2[2] - L1[2] * L2[1]
    Dx = -L1[3] * L2[2] - L1[2] * -L2[3]
    Dy = L1[1] * -L2[3] - -L1[3] * L2[1]
    if D != 0
        x = Dx / D
        y = Dy / D
        return [x,y]
    else
        error("no intersection")
    end
end

function get_circle_line_intersection(p1::Indexed_point, p2::Indexed_point, center::myPoint2d, radius::Float64)

    v = [p2.x-p1.x, p2.y-p1.y]

    # define t

    A = v[1]^2 + v[2]^2
    B = 2 * (v[1]*(p1.x-center.x)+v[2]*(p1.y-center.y))
    C = (p1.x-center.x)^2 + (p1.y-center.y)^2 - radius^2

    # above is the equation At^2 + Bt + C = 0;

    if B^2 - 4*A*C < 0

        return nothing, nothing
    end


    t1 = (-B - sqrt(B^2-4*A*C)) / (2 * A)
    t2 = (-B + sqrt(B^2-4*A*C)) / (2 * A)


    return t1, t2

end


function get_circle_line_intersection(p1::Indexed_point, p2::Indexed_point, center::Indexed_point, radius::Float64)

    v = [p2.x-p1.x, p2.y-p1.y]

    # define t

    A = v[1]^2 + v[2]^2
    B = 2 * (v[1]*(p1.x-center.x)+v[2]*(p1.y-center.y))
    C = (p1.x-center.x)^2 + (p1.y-center.y)^2 - radius^2

    # above is the equation At^2 + Bt + C = 0;

    if B^2 - 4*A*C < 0

        return nothing, nothing
    end


    t1 = (-B - sqrt(B^2-4*A*C)) / (2 * A)
    t2 = (-B + sqrt(B^2-4*A*C)) / (2 * A)


    return t1, t2

end

function get_t_from_xy(e::Indexed_edge, p::Indexed_point)
    return (p.x-e.s.x)/(e.e.x-e.s.x)
end

function get_xy_from_t(e::Indexed_edge, t::Float64)
    v = [e.e.x-e.s.x, e.e.y-e.s.y]
    p = [e.s.x, e.s.y]
    p = p+v*t
    return myPoint2d(p[1], p[2])
end

function get_xy_from_t(e::Indexed_edge, p::Indexed_point, t::Float64)
    v = [e.e.x-e.s.x, e.e.y-e.s.y]
    p = [e.s.x, e.s.y]
    p_new=p+v*t
    return Indexed_point(p.index, p_new[1], p_new[2])
end


function inward_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function signed_distance(e::myType.Indexed_edge, p::myType.Indexed_point)
    v = [p.x-e.s.x, p.y-e.s.y]
    n = inward_normal(e.s, e.e)
    return v[1]*n[1]+v[2]*n[2]
end


function signed_distance(p1, p2, p3)
    v = [p3.x-p1.x, p3.y-p1.y]
    n = inward_normal(p1, p2)
    return v[1]*n[1]+v[2]*n[2]
end

function compute_signed_angle(v1, v2)

	vn = [0,0,1]


	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return th

end

function inside_socket(p1::myType.Indexed_point, socket::myType.mySocket, epsilon = 0.00001)

	if p1.y >= socket.points[1].y
		return true
	end

	if point_on_socket(p1, socket, epsilon)
		return true
	end

	p = Luxor.Point(p1.x, p1.y)
	ps = Array{Luxor.Point,1}(undef, 0)
	for i = 1:length(socket.points)
		push!(ps, Luxor.Point(socket.points[i].x, socket.points[i].y))
	end

	return Luxor.isinside(p, ps)

end

function all_points_inside_socket(peg::myType.myPeg, socket::myType.mySocket, epsilon = 0.00001)
	l_p = length(peg.points)
	for i = 1:l_p
		if !inside_socket(peg.points[i], socket, epsilon)
			# println(i)
			return false
		end
	end

	return  true
end

function other_points_inside_socket(peg::myType.myPeg, socket::myType.mySocket, cm::myType.cm_id, epsilon = 0.00001)
	l_p = length(peg.points)
	p_ids = []
	for i = length(cm.cps)
		push!(p_ids, cm.cps[i].p_id)
	end

	for i = 1:l_p
		if peg.points[i].index in p_ids
			continue
		end


		if !inside_socket(peg.points[i], socket, epsilon)
			# println(i)
			return false
		end
	end

	return  true
end

function point_on_edge(p::myType.Indexed_point, e::myType.Indexed_edge, epsilon = 0.00001)
	Line = get_line_eq(e.s, e.e)
	if abs(Line[1]*p.x + Line[2]*p.y + Line[3]) <= epsilon && p.x >= min(e.s.x, e.e.x) && p.x <= max(e.s.x, e.e.x)
		return true
	else
		return false
	end
end

function point_on_socket(p::myType.Indexed_point, socket::myType.mySocket, epsilon = 0.00001)
	es = myType.mySocketEdges(socket)

	for i = 1:length(es.edges)
		if point_on_edge(p, es.edges[i], epsilon)
			return true
		end
	end

	if (p.y == socket.points[1].y && p.x <= socket.points[1].x)&&(p.y == socket.points[length(socket.points)].y && p.x >= socket.points[length(socket.points)].x)
		return true
	end


	return false
end




function will_rotate_outside(p1::myType.Indexed_point, socket::myType.mySocket, center::myType.myPoint2d, counter_clockwise = true)
	if !counter_clockwise
		p = Kinematics.rotate_p1_around_p2(p1, center, -0.01)

		return !inside_socket(p, socket)
	else
		p = Kinematics.rotate_p1_around_p2(p1, center, 0.01)

		return !inside_socket(p, socket)
	end
end

function two_circle_intersection(c1, c2, r1, r2)
	if r1 <= r2
		r1 += 0.000000000001
	else
		r2 += 0.000000000001
	end

	return Luxor.intersectioncirclecircle(Luxor.Point(c1.x, c1.y), r1, Luxor.Point(c2.x, c2.y), r2)
end

function no_flip(c1, c2, p_cc, peg)
	p1 = peg.points[c1.index]
	p2 = peg.points[c2.index]
	p3 = peg.points[p_cc.index]

	sign1 = sign(signed_distance(p1, p2, p3))
	sign2 = sign(signed_distance(c1, c2, p_cc))

	if sign1 == sign2
		return true
	else
		return false
	end

end

function get_one_point_from_two_points(c1, c2, r1, r2, p_id, peg)
	flag, ip1, ip2 = two_circle_intersection(c1, c2, r1, r2)
	if !flag
		# println("here: ", [c1, c2, r1, r2])
		return nothing
	end

	p1 = myType.Indexed_point(p_id, ip1.x, ip1.y)
	p2 = myType.Indexed_point(p_id, ip2.x, ip2.y)

	if no_flip(c1, c2, p1, peg)
		return p1
	elseif no_flip(c1, c2, p2, peg)
		return p2
	end

	return nothing

end

function compute_configuration(peg1, peg2)

	v1 = [0-peg.points[1].x, 0-peg.points[1].y]
	o1 = [0,0]
	o2 = [peg2.points[1]]

end

function test3(i)
	peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
	c1 = peg.points[1]
	c2 = peg.points[3]
	p = peg.points[6]
	r1 = sqrt((p.x - c1.x)^2 + (p.y - c1.y)^2)
	r2 = sqrt((p.x - c2.x)^2 + (p.y - c2.y)^2)
	# println(p)
	get_one_point_from_two_points(c1, c2, r1, r2, p.index, peg)

end

function test()
    p1 = [1, 1]
    p2 = [0, 0]
    p3 = [3, 0]
    p4 = [0, 3]
    line1 = get_line_eq(p1,p2)
    line2 = get_line_eq(p3,p4)
    p = get_two_line_intersection(line1, line2)

    p5 = Indexed_point(1, 0., 0.)
    p6 = Indexed_point(2, 1., 1.)
    r = 1.
    c = myPoint2d(0, 0)

    t1, t2 = get_circle_line_intersection(p5, p6, c, r)
end

function test2(i)
	peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)

	will_rotate_outside(peg.points[1], socket, myType.myPoint2d(peg.points[3].x, peg.points[3].y), false)

end


end
