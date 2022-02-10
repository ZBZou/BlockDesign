module ContactMode

using ..myType
using Combinatorics
using Compose
using Gadfly
using Colors
using GR

import ..Construct_Peg_Socket
import ..Constraints
import ..Kinematics
import ..ConcaveTest
import ..Draw
import PyPlot

function all_close(a, b)
	if abs(a - b) < 0.001
		return true
	end
	return false;

end


function is_contact_pair_valid(peg::myType.myPeg, socket::myType.mySocket, cp::myType.contact_pair)
	p_id = cp.c_p.index
	e_id = cp.c_e.index

    c_p = peg.points[p_id]
    es = mySocketEdges(socket)
    c_e = es.edges[e_id]

    np = length(peg.points)
	ne = length(socket.edges)

	result1 = false
	new_peg = nothing
	valid = []
	config_list = Array{myType.configSE2}(undef,0)
	th_temp = []
    for t = 0:0.01:1
		# println(t)
        pt = Constraints.get_xy_from_t(c_e, t)
        result, th = collision_check(pt, peg, socket, cp)

		if !result
			result1 = true
			for thi in th
				if !(thi in th_temp)
					push!(th_temp, thi)
					x = pt.x - c_p.x
					y = pt.y - c_p.y
					o = Kinematics.rotate_p1_around_p2([x,y],[pt.x, pt.y], thi)
					config = myType.configSE2(o[1],o[2], thi)
					push!(config_list, config)
				end
			end
			break
		end

		# if !result
		# 	result1 = true
		#
		# 	# x = pt.x - c_p.x
		# 	# y = pt.y - c_p.y
		# 	#
		# 	# config = myType.configSE2(x,y, 0)
		# 	# new_peg = Kinematics.get_peg_from_config(peg, config)
		#
		# 	if abs(th[1] - 0) <= abs(2*pi - th[2])
		# 		push!(valid, [abs(th[1] - 0), th[1], pt, th[2]])
		# 		# min_th = th[1]
		# 	else
		# 		# min_th = th[2]
		# 		push!(valid, [abs(2*pi - th[2]), th[2], pt, th[1]])
		# 	end
		#
		# 	# Draw.draw_sth(new_peg)
		# 	# new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
		#
		#
        #     # return true, new_peg
        # end

    end

	if !result1
		return result1, nothing
	end

	sort!(config_list, by = x -> x.th)
	# sort!(valid, by = x -> x[1])
	#
	# pt = valid[1][3]
	# min_th = valid[1][2]
	# max_th = valid[length(valid)][2]
	# pt2 = valid[length(valid)][3]
	#
	# x = pt.x - c_p.x
	# y = pt.y - c_p.y
	# x2 = pt2.x - c_p.x
	# y2 = pt2.y - c_p.y
	#
	# config = myType.configSE2(x,y, 0)
	# new_peg = Kinematics.get_peg_from_config(peg, config)
	# new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
	# o = Kinematics.rotate_p1_around_p2([x2,y2],[pt2.x, pt2.y], max_th)
	# new_config = myType.configSE2(o[1],o[2], max_th)
	# println([x,y],[o])
	#
    # return result1, new_peg, new_config

	return result1, config_list

end

function clockwise_angle(p1::Indexed_point, center::myPoint2d, e::Indexed_edge, cp::contact_pair)
	r = sqrt((p1.x - center.x)^2+(p1.y -center.y)^2)

	t1, t2 = Constraints.get_circle_line_intersection(e.s, e.e, center, r)


	if Constraints.point_on_edge(p1, e)
		if t1 >= 0 && t1 <= 1 && t2 >=0 && t2 <= 1
			return [0, 2*pi]
		else
			return [0]
		end
	end

	if (t1 == nothing) || ((t1 < 0 || t1 > 1) && (t2 < 0 || t2 > 1))
		return nothing
	end

	pt1 = Constraints.get_xy_from_t(e, t1)
	pt2 = Constraints.get_xy_from_t(e, t2)

	v1 = [p1.x - center.x, p1.y - center.y, 0]
	v2 = [pt1.x - center.x, pt1.y - center.y, 0]
	v3 = [pt2.x - center.x, pt2.y - center.y, 0]

	th1 = Constraints.compute_signed_angle(v1, v2)
	th2 = Constraints.compute_signed_angle(v1, v3)

	if th1 > 0
		th1 = th1 - 2*pi
	end

	if th2 > 0
		th2 = th2 - 2*pi
	end

	if ( 0<= t1 <=1) && (0<= t2 <=1)

		return [th1, th2]

	elseif ( 0<= t1 <=1) && (t2 < 0 || t2 > 1)

		return [th1]

	elseif ( 0<= t2 <=1) && (t1 < 0 || t1 > 1)
		return [th2]
	end

end

function counter_clockwise_angle(p1::Indexed_point, center::myPoint2d, e::Indexed_edge, cp::contact_pair)
	r = sqrt((p1.x - center.x)^2+(p1.y -center.y)^2)


	t1, t2 = Constraints.get_circle_line_intersection(e.s, e.e, center, r)



	if Constraints.point_on_edge(p1, e)
		if t1 >= 0 && t1 <= 1 && t2 >=0 && t2 <= 1
			return [0, 2*pi]
		else
			return [0]
		end
	end

	if (t1 == nothing) || ((t1 < 0 || t1 > 1) && (t2 < 0 || t2 > 1))
		return nothing
	end

	pt1 = Constraints.get_xy_from_t(e, t1)
	pt2 = Constraints.get_xy_from_t(e, t2)



	v1 = [p1.x - center.x, p1.y - center.y, 0]
	v2 = [pt1.x - center.x, pt1.y - center.y, 0]
	v3 = [pt2.x - center.x, pt2.y - center.y, 0]

	th1 = Constraints.compute_signed_angle(v1, v2)
	th2 = Constraints.compute_signed_angle(v1, v3)
	# println([e.index, pt2, th2])

	if th1 < 0
		th1 = 2*pi + th1
	end

	if th2 < 0
		th2 = 2*pi + th2
	end

	if ( 0<= t1 <=1) && (0<= t2 <=1)

		return [th1, th2]

	elseif ( 0<= t1 <=1) && (t2 < 0 || t2 > 1)

		return [th1]

	elseif ( 0<= t2 <=1) && (t1 < 0 || t1 > 1)
		return [th2]
	end
end

# function point_clockwise_rotation_limit(p1::Indexed_point, pt::myPoint2d, socket::mySocket, cp::contact_pair)
#     es = mySocketEdges(socket)
#     ne = length(es.edges)
#     th_least = []
#     th_most = []
# 	clockwise = true
#     for i = 1:ne
#         if !valid_contact(p1.index, i, cp)
#             continue
#         end
#
#         thi = clockwise_angle(p1, pt, es.edges[i], cp)
#         if thi == nothing
#             continue
#         end
#
# 		if !Constraints.inside_socket(p1, socket)
#             push!(th_least, thi)
#         elseif thi != 0
#             push!(th_most, thi)
# 		elseif thi == 0
# 			if Constraints.will_rotate_outside(p1, socket, pt, clockwise)
# 				push!(th_most, 0)
# 			end
# 		end
#
#     end
#
# 	if th_least == []
# 		th1 = nothing
# 	else
# 		sort!(th_least)
# 		th1 = th_least[1]
# 	end
#
# 	if th_most == []
# 		th2 = nothing
# 	else
# 		sort!(th_most)
# 		th2 = th_most(length(th_most)
# 	end
#
#     return th_least, th_most
# end
function point_clockwise_rotation_limit(p1::Indexed_point, pt::myPoint2d, socket::mySocket, cp::contact_pair)
    es = mySocketEdges(socket)
    ne = length(es.edges)
    # th_least = []
    # th_most = []
	th = []
	e_f = myType.Indexed_edge(0, socket.points[1], myType.Indexed_point(0, -10000000, socket.points[1].y))
	e_l = myType.Indexed_edge(ne+1, socket.points[length(socket.points)], myType.Indexed_point(ne+1, 10000000, socket.points[1].y))

	push!(es.edges, e_f)
	push!(es.edges, e_l)

    for i = 0:ne+1

		for e_i in es.edges
			if e_i.index == i
				thi =  clockwise_angle(p1, pt, e_i, cp)
			end
		end

		if thi == nothing
			continue
		end

		for th_j in thi
			push!(th, th_j)
		end
	end

	return th


		# if !Constraints.inside_socket(p1, socket)
        #     push!(th_least, thi)
        # elseif thi != 0
        #     push!(th_most, thi)
		# elseif thi == 0
		# 	if Constraints.will_rotate_outside(p1, socket, pt, clockwise)
		# 		push!(th_most, 0)
		# 	end
		# end
	#
    # end
	#
    # return th_least, th_most
end

function point_counter_clockwise_rotation_limit(p1::Indexed_point, pt::myPoint2d, socket::mySocket, cp::contact_pair)
	es = mySocketEdges(socket)
    ne = length(es.edges)
    # th_least = []
    # th_most = []
	th = []
	e_f = myType.Indexed_edge(0, socket.points[1], myType.Indexed_point(0, -10000000., socket.points[1].y))
	e_l = myType.Indexed_edge(ne+1, socket.points[length(socket.points)], myType.Indexed_point(ne+1, 10000000., socket.points[length(socket.points)].y))

	push!(es.edges, e_f)
	push!(es.edges, e_l)

    for i = 0:ne+1

		for e_i in es.edges
			if e_i.index == i
				thi =  counter_clockwise_angle(p1, pt, e_i, cp)

				if thi == nothing
					continue
				end

				for th_j in thi
					push!(th, th_j)
				end
			end
		end


	end

	return th


		# if !Constraints.inside_socket(p1, socket)
        #     push!(th_least, thi)
        # elseif thi != 0
        #     push!(th_most, thi)
		# elseif thi == 0
		# 	if Constraints.will_rotate_outside(p1, socket, pt, clockwise)
		# 		push!(th_most, 0)
		# 	end
		# end
	#
    # end
	#
    # return th_least, th_most
end

function valid_clockwise_rotation_angle(pt::myPoint2d, peg::myPeg, socket::mySocket, cp::contact_pair, epsilon = 0.001)
	np = length(peg.points)
	valid_th = []
	for i = 1:np
		valid_th_i = []

		if peg.points[i].index == cp.c_p.index
			continue
		end
		p_i = peg.points[i]
        th = point_clockwise_rotation_limit(p_i, pt, socket, cp)
		if !Constraints.inside_socket(p1, socket)
			if th == []
				return nothing
			end
			n_th = length(th)
			th_j = Array{float64}(undef, 2)
			for j = 1:n_th
				if isodd(j)
					th_j[1] = th[j]
				else
					th_j[2] = th[j]
					push!(valid_th_i, th_j)
				end
			end
		else
			if th == []
				push!(valid_th_i, [0, 2*pi])
			end

			n_th = length(th)
			th_j = Array{float64}(undef, 2)
			for j = 1:n_th
				if j == 1
					th_j[1] = 0
					th_j[2] = th[j]
					push!(valid_th_i, th_j)

				elseif j == n_th
					th_j[1] = th[j]
					th_j[2] = 2*pi
					push!(valid_th_i, th_j)

				else
					if iseven(j)
						th_j[1] = th[j]
					else
						th_j[2] = th[j]
						push!(valid_th_i, th_j)
					end
				end
			end


			push!(valid_th, valid_th_i)
		end


	end
	return valid_th

end

function valid_counter_clockwise_rotation_angle(pt::myPoint2d, peg::myPeg, socket::mySocket, cp::contact_pair, epsilon = 0.001)
	np = length(peg.points)
	valid_th = []
	for i = 1:np
		valid_th_i = []
		if peg.points[i].index == cp.c_p.index
			continue
		end
		p_i = peg.points[i]
        th = point_counter_clockwise_rotation_limit(p_i, pt, socket, cp)
		sort!(th)
		unique!(th)


		if Constraints.point_on_socket(p_i, socket)

			if !Constraints.will_rotate_outside(p_i, socket, pt)
				n_th = length(th)
				th_j = Array{Float64}(undef, 2)
				for j = 1:n_th
					if isodd(j)
						th_j[1] = th[j]
					else
						th_j[2] = th[j]
						push!(valid_th_i, th_j)
					end
				end
			else
				n_th = length(th)

				th_j = Array{Float64}(undef, 2)
				for j = 1:n_th
					if j == 1
						continue
					end

					if j == n_th
						th_j[1] = th[j]
						th_j[2] = 2*pi
						push!(valid_th_i, th_j)
						th_j = Array{Float64}(undef, 2)

					else
						if iseven(j)
							th_j[1] = th[j]
						else
							th_j[2] = th[j]
							push!(valid_th_i, th_j)
							th_j = Array{Float64}(undef, 2)
						end
					end


				end
			end


		else
			if !Constraints.inside_socket(p_i, socket)

				if th == []
					return nothing
				end
				n_th = length(th)
				th_j = Array{Float64}(undef, 2)
				for j = 1:n_th
					if isodd(j)
						th_j[1] = th[j]
					else
						th_j[2] = th[j]
						push!(valid_th_i, th_j)
					end
				end
			else
				if th == []
					push!(valid_th_i, [0, 2*pi])
					continue
				end

				n_th = length(th)

				th_j = Array{Float64}(undef, 2)
				for j = 1:n_th

					if j == 1
						th_j[1] = 0
						th_j[2] = th[j]

						push!(valid_th_i, th_j)

						th_j = Array{Float64}(undef, 2)



					elseif j == n_th

						th_j[1] = th[j]

						th_j[2] = 2*pi

						push!(valid_th_i, th_j)

						th_j = Array{Float64}(undef, 2)


					else
						if iseven(j)
							th_j[1] = th[j]
						else
							th_j[2] = th[j]
							push!(valid_th_i, th_j)
							th_j = Array{Float64}(undef, 2)
						end

					end

				end
			end


		end
		push!(valid_th, valid_th_i)
		# println([i, valid_th_i])

	end
	return valid_th


end

function collision_check(pt::myPoint2d, peg::myPeg, socket::mySocket, cp::contact_pair, epsilon = 0)
    np = length(peg.points)
	# println("here")

	config = myType.configSE2(pt.x - cp.c_p.x, pt.y - cp.c_p.y, 0)

	new_peg = Kinematics.get_peg_from_config(peg, config)
	#
	# Draw.draw_sth(new_peg)
	# Draw.draw_open_sth(socket, (0,0,1))



	th_cc = []
	l_list = []

	th_cc = valid_counter_clockwise_rotation_angle(pt, new_peg, socket, cp)

	# println("th_cc: ", th_cc)

	
	if th_cc == nothing
		return true, nothing
	end

	n_th_cc = length(th_cc)
	
	if n_th_cc == 0
		return true, nothing
	end

	count = 1

	for i = 1:n_th_cc
		if th_cc[i] == nothing
			return true, nothing
		end

		count *= length(th_cc[i])
		push!(l_list, length(th_cc[i]))
	end

	seq = find_combination(l_list)


	n_seq = length(seq)


	for i = 1:n_seq
		min_th = []
		max_th = []
		for j = 1:n_th_cc

			push!(min_th, correct_angle(th_cc[j][seq[i][j]][1]))
			push!(max_th, correct_angle(th_cc[j][seq[i][j]][2]))
		end
		sort!(min_th)
		sort!(max_th)
		# println(min_th)
		# println(max_th)
		
		# println("min_th: ", min_th)
		# println("max_th: ", max_th)

		if min_th[length(min_th)]+epsilon < max_th[1]-epsilon
			# println("here1")
			# println([min_th[length(min_th)]+epsilon, max_th[1]-epsilon])
			if min_th[length(min_th)]+epsilon <= pi/2 && max_th[1]-epsilon <= pi/2
				return false, [min_th[length(min_th)]+epsilon, max_th[1]-epsilon]
			elseif min_th[length(min_th)]+epsilon <= pi/2 && max_th[1]-epsilon > pi/2
				return false, [min_th[length(min_th)]+epsilon, pi/4]
			elseif max_th[1]-epsilon >= pi*3/2 && max_th[1]-epsilon >= pi*3/2
				return false, [min_th[length(min_th)]+epsilon, max_th[1]-epsilon]
			elseif max_th[1]-epsilon >= pi*3/2 && max_th[1]-epsilon < pi*3/2
				return false, [pi*3/2, max_th[1]-epsilon]
			end
			# println("here2")
		end
	end


	return true, nothing

	#
    #     pi = peg.points[i]
    #     th_least, th_most = point_rotation_limit(pi, pt, socket, cp)
    #     # if i < cp.c_p.index
    #     #     for thi in th
    #     #         if thi > 0
    #     #             push!(th_right_least, thi)
    #     #         else
    #     #             push!(th_left_most, thi)
    #     #         end
    #     #     end
    #     # else
    #     #     for thi in th
    #     #         if thi > 0
    #     #             push!(th_right_most, thi)
    #     #         else
    #     #             push!(th_left_least, thi)
    #     #         end
    #     #     end
    #     # end
	#
    #     for thi in th_least
    #         if thi > 0
    #             push!(th_right_least, thi)
    #         else
    #             push!(th_left_least, thi)
    #         end
    #     end
	#
	#
    #     for thj in th_most
    #         if thj > 0
    #             push!(th_right_most, thj)
    #         else
    #             push!(th_left_most, thj)
    #         end
    #     end
	#
    # end
	#
    # sort!(th_right_least)
    # sort!(th_right_most)
    # sort!(th_left_least)
    # sort!(th_left_most)
    # println([th_right_least, th_right_most])
    # println([th_left_least, th_left_most])
	#
    # if th_right_least !=[] && th_right_most !=[] && (th_right_least[length(th_right_least)] > th_right_most[1])
    #     result1 = true
    # else
    #     result1 = false
    # end
	#
    # if th_left_least !=[] && th_left_most !=[] && (th_left_least[1] < th_left_most[length(th_left_most)])
    #     result2 = true
    # else
    #     result2 = false
    # end
	#
    # return result1 || result2
	#
	#
	#
	#
    # #     if pt.x <= peg.points[cp.c_p.index].x
    # #         if i < cp.c_p.index
    # #             push!(th1, thi[length(thi)])
    # #         else
    # #             push!(th2, thi[1])
    # #         end
    # #     else
    # #         if i > cp.c_p.index
    # #             push!(th1, thi[length(thi)])
    # #         else
    # #             push!(th2, thi[1])
    # #         end
    # #     end
    # # end
    # # println([th1, th2])
    # # sort!(th1)
    # # sort!(th2)
    # #
    # # if th2[1] - th1[length(th1)] > epsilon
    # #     return false
    # # else
    # #     return true
    # # end

end

function correct_angle(ang)
	while ang >= 2*pi || ang <= -2*pi
		if ang >= 2*pi
			ang = ang - 2*pi
		else
			ang += 2*pi
		end
	end
	return ang
end

function find_all_pairs(peg::myType.myPeg, socket::myType.mySocket, contacts::Array{contact})
	# for this, should select and test all the possible pairs

	np = length(peg.points)
	ne = length(socket.edges)
	es = myType.mySocketEdges(socket)
    n = ne/2
	cps = []
	pegs = []

	c_lists = []

	c_ids = get_initial_p_e_contact_id(peg, contacts)

    if ConcaveTest.is_concave(peg)
        p_concave = ConcaveTest.concave_points(peg)
    else
        p_concave = []
    end

    # count = 0
	for i = 1:np

        if i in p_concave
            # count += 1
            continue
        end
        # sq = i - count

		for j = 1:ne

			 if !valid_contact(i, j, socket, c_ids::Array)
				 # println([i, j])
				 continue
			 end
            # if (contacts[sq].e_id <= n && j> contacts[sq].e_id) || (contacts[sq].e_id >= n+1 && j < contacts[sq].e_id)
            #     continue
            # end



			# try each possible combination
			# test pair

			cp = myType.contact_pair((i-1)*ne + j, peg.points[i], es.edges[j])
			result, config_list = is_contact_pair_valid(peg, socket, cp)
			if !result
				continue
			end

			push!(cps, myType.contact_id(i, j))
			push!(c_lists, config_list)

		end
	end

	return cps, c_lists

	# println(ps)

end

function find_all_pairs(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{contact_id})
	# for this, should select and test all the possible pairs

	np = length(peg.points)
	ne = length(socket.edges)
	es = myType.mySocketEdges(socket)
    n = ne/2
	cps = []
	pegs = []
	c_lists = []

    if ConcaveTest.is_concave(peg)
        p_concave = ConcaveTest.concave_points(peg)
    else
        p_concave = []
    end

    # count = 0
	for i = 1:np

        if i in p_concave
            # count += 1
            continue
        end
        # sq = i - count

		for j = 1:ne

			 if !valid_contact(i, j, socket, c_ids::Array)
				 # println([i, j])
				 continue
			 end
            # if (contacts[sq].e_id <= n && j> contacts[sq].e_id) || (contacts[sq].e_id >= n+1 && j < contacts[sq].e_id)
            #     continue
            # end


			# try each possible combination
			# test pair

			cp = myType.contact_pair((i-1)*ne + j, peg.points[i], es.edges[j])
			result, config_list = is_contact_pair_valid(peg, socket, cp)
			if !result
				continue
			end

			push!(cps, myType.contact_id(i, j))
			push!(c_lists, config_list)

		end
	end
	# println(cps)
	return cps, c_lists

	# println(ps)

end

function find_combination(l_list::Array)
	n = length(l_list)
	count = 1
	c_list = []
	for i = 1:n
		count *= l_list[i]
		push!(c_list, count)
	end

	result = []

	for i = 1:count
		temp = Array{Int}(undef, n)
		m = i-1
		for j = n:-1:1
			if j == 1
				temp[j] =1+m
				continue
			end

			temp[j] = 1 + Int(floor(m/c_list[j-1]))
			m = m%c_list[j-1]
		end

		push!(result, temp)
	end

	return  result


end

function test3()
	find_combination([3,2,2])

end

# function valid_contact(p_id::Int, e_id::Int, cp::contact_pair)
#     if (p_id < cp.c_p.index && e_id <= cp.c_e.index)||(p_id > cp.c_p.index && e_id >= cp.c_e.index)
#         return true
#     else
#         return false
#     end
#
# end

function get_initial_p_e_contact_id(peg::myType.myPeg, cs::Array{contact})
	np = length(peg.points)

	c_ids = []

	if ConcaveTest.is_concave(peg)
		p_concave = ConcaveTest.concave_points(peg)
	else
		p_concave = []
	end

	count = 0
	for i = 1:np

		if i in p_concave
			count += 1
			continue
		end

		p_id = i

		sq = i - count

		e_id = cs[sq].e_id

		push!(c_ids, myType.contact_id(p_id, e_id))

	end

	return c_ids


end

function valid_contact(p_id::Int, e_id::Int, socket::myType.mySocket, c_ids::Array)
	l_cids = 0
	temp_c_ids = []
	for i = 1:length(c_ids)
		push!(temp_c_ids, c_ids[i].p_id)
	end
	unique!(temp_c_ids)
	l_cids = length(temp_c_ids)

	es = myType.mySocketEdges(socket)
	ne = length(socket.points) - 1
	# 0 mid, 1 right, -1 left
	p_label = 0
	e_label = 0
	if abs(p_id - 1) == abs(p_id - l_cids)
		p_label = 0
	elseif abs(p_id - 1) < abs(p_id - l_cids)
		p_label = -1
	else abs(p_id - 1) > abs(p_id - l_cids)
		p_label = 1
	end

	if abs(e_id - 1) == abs(e_id - ne)
		e_label = 0
	elseif abs(e_id - 1) < abs(e_id - ne)
		e_label = -1
	else abs(e_id - 1) > abs(e_id - ne)
		e_label = 1
	end


	for i = 1:length(c_ids)
		if p_id == c_ids[i].p_id
			id = c_ids[i].e_id

		else
			continue
		end

		if e_id == id
			# println([e_id, id, c_ids])
			return true
		end
		# println("es edges: ", length(es.edges), ", id: ", id)
		if id > length(es.edges)
			return false
		end
		e1 = es.edges[id]

		height1 = (e1.s.y + e1.e.y)/2

		e2 = es.edges[e_id]

		height2 = (e2.s.y + e2.e.y)/2

		if height2 > height1 && e_label * p_label >= 0
			return true
		end

	end

	return false

end

function is_pp_contact(cm::myType.cm_id)
	n = length(cm.cps)

	p_ids = []
	for i = 1:n
		push!(p_ids, cm.cps[i].p_id)
	end

	uni_id = unique!(p_ids)

	if length(uni_id) < n
		return true
	else
		return false
	end
end

function is_pp_contact(cm::myType.contact_mode)
	n = length(cm.cps)

	p_ids = []
	for i = 1:n
		push!(p_ids, cm.cps[i].c_p.index)
	end

	uni_id = unique!(p_ids)

	if length(uni_id) < n
		return true
	else
		return false
	end
end

function get_pp_id(cm::myType.cm_id, socket::myType.mySocket)
	n = length(cm.cps)
	id = []
	for i = 1:n
		if i == 1
			continue
		end

		last = i-1

		if cm.cps[i].p_id == cm.cps[last].p_id
			push!(id, [cm.cps[i].p_id, socket.edges[cm.cps[i].e_id].s])

		end

	end

	return id

end

function test_pp(i)
	peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
	cs = Construct_Peg_Socket.get_contacts(i)
	cps, c_list = find_all_pairs(peg, socket, cs)

	cm = myType.cm_id(cps)
	is_pp_contact(cm)
	get_pp_id(cm, socket)

end

function cm_with_sigle_cp(peg::myType.myPeg, socket::myType.mySocket, cm::myType.cm_id, epsilon=0.005)
	es = mySocketEdges(socket)

	p_id = cm.cps[1].p_id
	e_id = cm.cps[1].e_id

	c_p = peg.points[p_id]
	es = mySocketEdges(socket)
	c_e = es.edges[e_id]

	cp = myType.contact_pair(1, c_p, c_e)

	np = length(peg.points)
	ne = length(socket.edges)

	result1 = false

	valid = []
	th_temp = []
	config_list = Array{myType.configSE2}(undef,0)
	for t = 0.01:0.01:0.99
		# println(t)
        pt = Constraints.get_xy_from_t(c_e, t)
        result, th = collision_check(pt, peg, socket, cp, epsilon)

		if !result
			result1 = true
			for thi in th
				if !(thi in th_temp)
					push!(th_temp, thi)
					x = pt.x - c_p.x
					y = pt.y - c_p.y
					o = Kinematics.rotate_p1_around_p2([x,y],[pt.x, pt.y], thi)
					config = myType.configSE2(o[1],o[2], thi)
					push!(config_list, config)
				end
			end
		end

	end

	if !result1
		return result1, nothing
	end

	sort!(config_list, by = x -> x.th)


	return result1, config_list

	# 	if !result
	# 		result1 = true
	#
	# 		# x = pt.x - c_p.x
	# 		# y = pt.y - c_p.y
	# 		#
	# 		# config = myType.configSE2(x,y, 0)
	# 		# new_peg = Kinematics.get_peg_from_config(peg, config)
	#
	# 		if abs(th[1] - 0) <= abs(2*pi - th[2])
	# 			push!(valid, [abs(th[1] - 0), th[1], pt])
	# 			# min_th = th[1]
	# 		else
	# 			# min_th = th[2]
	# 			push!(valid, [abs(2*pi - th[2]), th[2], pt])
	# 		end
	#
	# 		# Draw.draw_sth(new_peg)
	# 		# new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
	#
	#
    #         # return true, new_peg
    #     end
	#
    # end
	#
	# if !result1
	# 	return result1, nothing
	# end
	#
	# sort!(valid, by = x -> x[1])
	#
	# pt = valid[1][3]
	# min_th = valid[1][2]
	#
	#
	# x = pt.x - c_p.x
	# y = pt.y - c_p.y
	#
	# config = myType.configSE2(x,y, 0)
	# new_peg = Kinematics.get_peg_from_config(peg, config)
	# new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
	#
    # return result1, new_peg
	#
	# # for t = 0.01:0.01:0.99
	# # 	# println(t)
	# # 	pt = Constraints.get_xy_from_t(c_e, t)
	# # 	result, th = collision_check(pt, peg, socket, cp, epsilon)
	# #
	# #
	# #
	# # 	if !result
	# #
	# # 		x = pt.x - c_p.x
	# # 		y = pt.y - c_p.y
	# # 		# mid_th = (th[1]+th[2])/2
	# # 		config = myType.configSE2(x,y, 0)
	# # 		new_peg = Kinematics.get_peg_from_config(peg, config)
	# #
	# # 		if abs(th[1] - 0) <= abs(2*pi - th[2])
	# # 			min_th = th[1]
	# # 		else
	# # 			min_th = th[2]
	# # 		end
	# # 		# println(min_th)
	# # 		# Draw.draw_sth(new_peg)
	# # 		new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
	# #
	# #
	# # 		return true, new_peg
	# # 	end
	# #
	# # end
	# # return false, nothing
end

function distance(p1, p2)
	return sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)
end

function distance_dic(peg::myType.myPeg)
	n = length(peg.points)
	dics = Dict()
	for i = 1:n
		dic_i = Dict()
		p_i = peg.points[i]
		for j = 1:n
			p_j = peg.points[j]
			dis = distance(p_i, p_j)
			get!(dic_i, j, dis)
		end
		get!(dics, i, dic_i)

	end
	return dics
end

function get_peg_from_two_points(peg, dis_dic, current)
	l_cur = length(current)
	l_p = length(peg.points)

	new_peg = deepcopy(peg)

	if l_cur < 2
		return nothing
	end
	p1 = []
	p2 = []
	count = 0
	for i = 1:l_p
		if haskey(current, i)
			count += 1

			if count == 1
				p = get(current, i, 1)
				p1 = myType.Indexed_point(i, p.x, p.y)
			elseif count == 2
				p = get(current, i, 1)
				p2 = myType.Indexed_point(i, p.x, p.y)
				break
			end
		end
	end

	for i = 1:l_p
		p_i = new_peg.points[i]
		if p_i.index == p1.index
			new_peg.points[i] = p1
			continue
		elseif p_i.index == p2.index
			new_peg.points[i] = p2
			continue
		end

		r1 = get(get(dis_dic, i, 1), p1.index, 1)
		r2 = get(get(dis_dic, i, 1), p2.index, 1)

		p_temp = Constraints.get_one_point_from_two_points(p1, p2, r1, r2, i, peg)

		if p_temp == nothing
			return nothing
		else
			new_peg.points[i] = p_temp
		end
	end
	return new_peg


end

function get_peg_rotation_th(peg, new_peg)
	p1 = peg.points[1]
	p2 = peg.points[2]
	np1 = new_peg.points[1]
	np2 = new_peg.points[2]

	v1 = [p2.x - p1.x, p2.y - p1.y, 0]

	v2 = [np2.x - np1.x, np2.y - np1.y, 0]

	th = Constraints.compute_signed_angle(v1, v2)
	return th
end

function re_sort_cm(cm::myType.cm_id)
	cps = []
	dis = []
	n = length(cm.cps)
	if n<=2
		return cm
	end

	for i = 1:n
		for j = 1:n
			if i == j || [i, j] in cps || [j,i] in cps
				continue
			end
			temp = abs(cm.cps[i].p_id - cm.cps[j].p_id)
			push!(dis, [i, j, temp])
			push!(cps, [i, j])
		end
	end

	sort!(dis, by=x->x[3])
	# println(dis)

	temp_cp=[]
	temp_cm=[]
	for i = 1:length(dis)
		for j = 1:2
			push!(temp_cp, dis[i][j])
		end
	end
	unique!(temp_cp)
	for i = 1:length(temp_cp)
		push!(temp_cm, cm.cps[temp_cp[i]])
	end

	return myType.cm_id(temp_cm)

end

function cm_with_mlti_cp(peg::myType.myPeg, socket::myType.mySocket, cm::myType.cm_id, dis_epsilon = 0.001)
	l_cps = length(cm.cps)
	index = 1
	cm = re_sort_cm(cm)
	# println(cm)
	c_e1 = cm.cps[index].e_id
	c_p1 = cm.cps[index].p_id
	result = false
	result1 = false
	dis_dic = distance_dic(peg)
	es = myType.mySocketEdges(socket)
	new_peg = nothing
	config_list = Array{myType.configSE2}(undef,0)
	valid = []
	p_dis = 0
	delta_t = 0.001
	sign = 1
	t = 0.001
	while t <=0.99 && delta_t>0.00001
		pt = Constraints.get_xy_from_t(es.edges[c_e1], t)
		current = Dict()
		get!(current, cm.cps[index].p_id, pt)
		result1, current = multi_circle_intersection(peg, socket, current, cm, index+1, dis_dic, result1)
		if !result1
			if current != nothing
				# println("here: ", current)
				if current * p_dis >= 0
					p_dis = current
					t = t+sign*delta_t
					continue
				else
					# println("here1: ", [delta_t, p_dis])
					sign = -sign
					delta_t = delta_t/2
					t = t+sign*delta_t
					p_dis = current
					# if delta_t<=0.000001
					# 	# println("here1")
					# 	@goto valid
					# end
					continue
				end
			end
		else
			# println(t)
			# println(current)
			new_peg = get_peg_from_two_points(peg, dis_dic, current)
			# Draw.draw_all_peg([new_peg], socket)
			# println([Constraints.other_points_inside_socket(new_peg, socket, cm, dis_epsilon), other_points_do_not_contact(current, new_peg, socket)])

			if Constraints.other_points_inside_socket(new_peg, socket, cm, dis_epsilon) && other_points_do_not_contact(current, new_peg, socket)
					# println("here2")
				result = true
				th = get_peg_rotation_th(peg, new_peg)
				x = pt.x - peg.points[c_p1].x
				y = pt.y - peg.points[c_p1].y
				o = Kinematics.rotate_p1_around_p2([x,y],[pt.x, pt.y], th)
				config = myType.configSE2(o[1],o[2],th)
				# push!(valid, [abs(th), new_peg])
				push!(config_list, config)
				break
			end
			delta_t = 0.001
			sign = 1

		end
		p_dis = 0
		t = t+sign*delta_t
	end

	if !result
		return result, new_peg
	end

	# sort!(valid, by = x->x[1])
	sort!(config_list, by = x->x.th)
	return result, config_list

	# return result, valid[1][2]

end

function other_points_do_not_contact(current::Dict, peg, socket)
	l_p = length(peg.points)

	p_concave = []

	if ConcaveTest.is_concave(peg)
		p_concave = ConcaveTest.concave_points(peg)
	else
		p_concave = []
	end
	# println(p_concave)

	for i = 1:l_p
		p_i = peg.points[i]
		if haskey(current, p_i.index) || p_i.index in p_concave
			# println(i)
			continue
		end

		if Constraints.point_on_socket(p_i, socket)
			# println(i)
			return false
		end
	end

	return true
end

function other_points_do_not_contact(cp_ids::Array, peg, socket)
	l_p = length(peg.points)

	p_concave = []

	if ConcaveTest.is_concave(peg)
		p_concave = ConcaveTest.concave_points(peg)
	else
		p_concave = []
	end
	# println(p_concave)

	for i = 1:l_p
		p_i = peg.points[i]
		if p_i.index in cp_ids || p_i.index in p_concave
			# println(i)
			continue
		end

		if Constraints.point_on_socket(p_i, socket)
			# println(i)
			return false
		end
	end

	return true
end

function multi_circle_intersection(peg, socket, current, cm, index, dis_dic, result, th_epsilon=0.001, dis_epsilon=0.005)
	es = myType.mySocketEdges(socket)


	if index > length(cm.cps)
		# println("here2")
		return true, current
	end

	cp = cm.cps[index]
	last_cp = cm.cps[index-1]

	p_id = cp.p_id
	e_id = cp.e_id
	last_pid = last_cp.p_id
	e = es.edges[e_id]

	if index == 2
		r = get(get(dis_dic, p_id, 1), last_pid, 1)
		center = get(current, last_pid, 1)
		t1, t2 = Constraints.get_circle_line_intersection(e.s, e.e, center, r)
		if t1 == nothing && t2 == nothing
			return false, nothing
		end

		if (t2>=th_epsilon && t2<= 1-th_epsilon) && (t1>=th_epsilon && t1<= 1-th_epsilon)
			pt = Constraints.get_xy_from_t(e, t1)
			cur = deepcopy(current)
			get!(cur, p_id, pt)
			result, cur = multi_circle_intersection(peg, socket, cur, cm, index+1, dis_dic, result)
			# println([t1,",", t2,",", result ])
			if !result
				pt = Constraints.get_xy_from_t(e, t2)
				cur = deepcopy(current)
				get!(cur, p_id, pt)
				result, cur = multi_circle_intersection(peg, socket, cur, cm, index+1, dis_dic, result)
			end
		elseif (t2>=th_epsilon && t2<= 1-th_epsilon) && (t1 < th_epsilon || t1 > th_epsilon)
			pt = Constraints.get_xy_from_t(e, t2)
			# println(pt)
			# println(Constraints.inside_socket(myType.Indexed_point(1, pt.x, pt.y), socket))
			#
			# println(Constraints.point_on_edge(myType.Indexed_point(1, pt.x, pt.y), e))
			cur = deepcopy(current)
			get!(cur, p_id, pt)
			result, cur = multi_circle_intersection(peg, socket, cur, cm, index+1, dis_dic, result)
		elseif (t1>=th_epsilon && t1<= 1-th_epsilon) && (t2 < th_epsilon || t2 > th_epsilon)
			pt = Constraints.get_xy_from_t(e, t1)
			# println(pt)
			cur = deepcopy(current)
			get!(cur, p_id, pt)
			result, cur = multi_circle_intersection(peg, socket, cur, cm, index+1, dis_dic, result)
		else
			return false, nothing
		end

	else
		l_p = length(dis_dic)
		p1 = []
		p2 = []
		count = 0
		for i = 1:l_p
			if haskey(current, i)
				count += 1

				if count == 1
					p = get(current, i, 1)
					p1 = myType.Indexed_point(i, p.x, p.y)
				elseif count == 2
					p = get(current, i, 1)
					p2 = myType.Indexed_point(i, p.x, p.y)
					break
				end
			end
		end

		r1 = get(get(dis_dic, p_id, 1), p1.index, 1)
		r2 = get(get(dis_dic, p_id, 1), p2.index, 1)

		new_peg = get_peg_from_two_points(peg, dis_dic, current)
		# Draw.draw_all_peg([new_peg], socket)

		p_temp = Constraints.get_one_point_from_two_points(p1, p2, r1, r2, p_id, peg)


		# println(e)
		if p_temp == nothing
			return false, nothing
		end


		# println(Constraints.point_on_edge(p_temp, e))
		if Constraints.point_on_edge(p_temp, e, dis_epsilon)
			# println("here", ",", index)
			# println("here")
			cur = deepcopy(current)
			get!(cur, p_id, p_temp)
			result, cur = multi_circle_intersection(peg, socket, cur, cm, index+1, dis_dic, result)
		elseif [p_temp.x - e.s.x, p_temp.y-e.s.y]'*[e.e.x - e.s.x, e.e.y-e.s.y]>0

			return false, Constraints.signed_distance(e, p_temp)
		else
			return false, nothing
		end

	end






	return result, cur

end

function cm_with_pp_cp(peg, socket, cm, th_epsilon=0.001, dis_epsilon=0.005)
	n = length(cm.cps)

	p_ids = []
	for i = 1:n
		push!(p_ids, cm.cps[i].p_id)
	end

	sort!(p_ids)

	pp_ids =[]

	for i = 1:n-1
		next = i+1
		if next == n+1
			next = 1
		end
		if p_ids[i] == p_ids[next]
			push!(pp_ids, p_ids[i])
		end
	end

	if length(pp_ids) == 0
		error("no pp")
	end

	if length(pp_ids) == 1
		result, c_list = cm_with_one_pp(peg, socket, cm, pp_ids[1], p_ids, th_epsilon)
	else
		result, c_list = cm_with_multi_pp(peg, socket, cm, pp_ids, p_ids, dis_epsilon)
	end

	return result, c_list
end


function cm_with_one_pp(peg, socket, cm, pp_id, p_ids, th_epsilon=0.001, dis_epsilon=0.005)
	e_ids = []
	n = length(cm.cps)
	for i = 1:n
		if pp_id == cm.cps[i].p_id
			push!(e_ids, cm.cps[i].e_id)
		end
	end
	sort!(e_ids)
	p_ids_temp = deepcopy(p_ids)
	unique!(p_ids_temp)
	# println([p_ids_temp, pp_id])
	if length(p_ids_temp) == 1 && pp_id in p_ids_temp

		c_p = peg.points[pp_id]
		es = mySocketEdges(socket)
		c_e = es.edges[e_ids[1]]

		cp = myType.contact_pair(1, c_p, c_e)

		result1 = false

		valid = []
		config_list = Array{myType.configSE2}(undef,0)
		th_temp = []
		pt = myPoint2d(c_e.e.x, c_e.e.y)
		result, th = collision_check(pt, peg, socket, cp, th_epsilon)
		if !result
			result1 = true
			for thi in th
				if !(thi in th_temp)
					push!(th_temp, thi)
					x = pt.x - c_p.x
					y = pt.y - c_p.y
					o = Kinematics.rotate_p1_around_p2([x,y],[pt.x, pt.y], thi)
					config = myType.configSE2(o[1],o[2], thi)
					push!(config_list, config)
				end
			end
		end

		sort!(config_list, by = x -> x.th)
		return result1, config_list

		# if !result
		# 	result1 = true
		#
		# 	# x = pt.x - c_p.x
		# 	# y = pt.y - c_p.y
		# 	#
		# 	# config = myType.configSE2(x,y, 0)
		# 	# new_peg = Kinematics.get_peg_from_config(peg, config)
		#
		# 	if abs(th[1] - 0) <= abs(2*pi - th[2])
		# 		push!(valid, [abs(th[1] - 0), th[1], pt])
		# 		# min_th = th[1]
		# 	else
		# 		# min_th = th[2]
		# 		push!(valid, [abs(2*pi - th[2]), th[2], pt])
		# 	end
		# end
		#
		# 	# Draw.draw_sth(new_peg)
		# 	# new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
		#
		#
		# 	# return true, new_peg
		#
		#
		# if !result1
		# 	return result1, nothing
		# end
		#
		# sort!(valid, by = x -> x[1])
		#
		# pt = valid[1][3]
		# min_th = valid[1][2]
		#
		#
		# x = pt.x - c_p.x
		# y = pt.y - c_p.y
		#
		# config = myType.configSE2(x,y, 0)
		# new_peg = Kinematics.get_peg_from_config(peg, config)
		# new_peg = Kinematics.rotate_peg_around_p(new_peg, pt, min_th)
		#
		# return result1, new_peg
	else
		config_list = Array{myType.configSE2}(undef,0)
		c_p = peg.points[pp_id]
		es = mySocketEdges(socket)
		c_e = es.edges[e_ids[1]]

		cp = myType.contact_pair(1, c_p, c_e)

		result1 = false

		pt = myPoint2d(c_e.e.x, c_e.e.y)

		index = 1

		result1 = false

		dis_dic = distance_dic(peg)

		new_peg = nothing

		# cp_invilid = myType.contact_pair(2, c_p, es.edges[e_ids[2]])
		cp_vilid = myType.contact_id(pp_id, e_ids[1])

		cm_temp = deepcopy(cm)

		deleteat!(cm_temp.cps, findall(x-> x.p_id==pp_id, cm_temp.cps))
		insert!(cm_temp.cps, 1, cp_vilid)



		current = Dict()
		get!(current, pp_id, pt)
		result1, current = multi_circle_intersection(peg, socket, current, cm_temp, index+1, dis_dic, result1)

		if !result1
			return false, nothing
		else
			# println(t)
			# println(current)
			new_peg = get_peg_from_two_points(peg, dis_dic, current)

			# Draw.draw_all_peg([new_peg], socket)
			# println([Constraints.all_points_inside_socket(new_peg, socket, dis_epsilon), other_points_do_not_contact(current, new_peg, socket)])

			if Constraints.all_points_inside_socket(new_peg, socket, dis_epsilon) && other_points_do_not_contact(current, new_peg, socket)
				th = get_peg_rotation_th(peg, new_peg)

				x = pt.x - peg.points[pp_id].x
				y = pt.y - peg.points[pp_id].y
				o = Kinematics.rotate_p1_around_p2([x,y],[pt.x,pt.y],th)
				config = myType.configSE2(o[1],o[2],th)
				push!(config_list, config)
				return true, config_list

			else
				return false, nothing
			end
		end

	end

end

function cm_with_multi_pp(peg, socket, cm, pp_ids, p_ids, dis_epsilon = 0.005)


	es = myType.mySocketEdges(socket)
	dis_dic = distance_dic(peg)
	displacement = Dict()
	config_list = Array{myType.configSE2}(undef,0)
	for i = 1:2
		e_ids_i = []
		for j = 1:length(cm.cps)
			if pp_ids[i] == cm.cps[j].p_id
				push!(e_ids_i, cm.cps[i].e_id)
			end
		end
		sort!(e_ids_i)
		get(displacement, pp_ids[i], es.edges[e_ids_i[1]].e)
	end

	new_peg = get_peg_from_two_points(peg, dis_dic, displacement)
	es = myType.mySocketEdges(socket)

	if new_peg == nothing
		return false, nothing
	end

	l_cps = length(cm.cps)

	for i = 1:l_cps
		if cm.cps[i].p_id == pp_ids[1] || cm.cps[i].p_id == pp_ids[2]
			continue
		end

		p_i = new_peg.points[cm.cps[i].p_id]
		e_i = es.edges[cm.cps[i].e_id]

		if !Constraints.point_on_edge(p_i, e_i, dis_epsilon)
			return false, nothing
		end

	end

	if Constraints.all_points_inside_socket(new_peg, socket, dis_epsilon) && other_points_do_not_contact(p_ids, new_peg, socket)
		th = get_peg_rotation_th(peg, new_peg)
		x = new_peg.point[pp_ids[1]].x - peg.points[pp_ids[1]].x
		y = new_peg.point[pp_ids[1]].y - peg.points[pp_ids[1]].y
		o = Kinematics.rotate_p1_around_p2([x,y],[new_peg.point[pp_ids[1]].x,new_peg.point[pp_ids[1]].y],th)


		config = myType.configSE2(o[1],o[2],th)
		push!(config_list, config)
		return true, config_list
	else
		return false, nothing
	end

end

function is_cm_valid(peg::myType.myPeg, socket::myType.mySocket, cm::myType.cm_id)

	if length(cm.cps) == 1

		return  cm_with_sigle_cp(peg, socket, cm)
	else
		if !is_pp_contact(cm)
			return cm_with_mlti_cp(peg, socket, cm)
		else
			return cm_with_pp_cp(peg, socket, cm)
		end
	end

end


function possible_contact_modes(peg::myType.myPeg, socket::myType.mySocket, cs::Array{contact})
	l_p = length(peg.points)
	l_e = length(socket.edges)
	c_ids = get_initial_p_e_contact_id(peg, cs)

	cps,c_list = find_all_pairs(peg, socket, cs)

	l_cps = length(cps)

	# println(cps)

	cms = []


	for i = 1:l_cps

		cm_i = pick_k(i, cps)

		if length(cm_i) != 0
			push!(cms, cm_i)
		end


	end

	return cms


end


function possible_contact_modes(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{contact_id})
	l_p = length(peg.points)
	l_e = length(socket.edges)

	cps,c_list = find_all_pairs(peg, socket, c_ids)

	l_cps = length(cps)

	# println(cps)

	cms = []


	for i = 1:l_cps

		cm_i = pick_k(i, cps)

		if length(cm_i) != 0
			push!(cms, cm_i)
		end


	end

	return cms


end


function possible_three_pairs_cms(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{contact_id})
	l_p = length(peg.points)
	l_e = length(socket.edges)

	cps,c_list = find_all_pairs(peg, socket, c_ids)

	l_cps = length(cps)

	# println(cps)

	cms = []

	for i = 3:l_cps

		cm_i = pick_k(i, cps)

		if length(cm_i) != 0
			push!(cms, cm_i)
		end


	end

	return cms
end

function validate_new_pair(cps, current, newcp)


	if length(current) == 1

		cpi = cps[current[1]].p_id
		cei = cps[current[1]].e_id

		cpj = cps[newcp].p_id
		cej = cps[newcp].e_id

		if cpi == cpj && abs(cei - cej) != 1

			return false
		end

		# if cei == cej && abs(cpi - cpj) != 1
		#
		# 	return false
		# end

		return true
	end


	c_p = cps[newcp].p_id
	c_e = cps[newcp].e_id

	appearance = Dict{Int, Int}()

	for i = 1:length(current)
		cpi = cps[current[i]].p_id
		cei = cps[current[i]].e_id

		if haskey(appearance, cpi)
			appearance[cpi] += 1
		else
			appearance[cpi] = 1
		end

		if cpi == c_p && abs(cei - c_e) != 1
			return false
		end

		if cpi == c_p && appearance[cpi] > 1
			return false
		end

		# if cei == c_e && abs(cpi - c_p) != 1
		#
		# 	return false
		# end

	end

	return true

end


function pick_k_from_array(k::Int, indicies::Array{Int}, result::Array{Int}, total, cps)

	# println(indicies)
	# println(k)
	if k == 0
		push!(total, result)
	end

	for i = 1:length(indicies)

		if indicies[i] < result[length(result)]
			continue
		end

		if ! validate_new_pair(cps, result, indicies[i])

			continue
		end

		nids = deepcopy(indicies)
		nrs = deepcopy(result)

		push!(nrs, nids[i])
		# println([indicies[i], nrs])

		# println("before",nids, ",",nids[i])
		filter!(x->x≠indicies[i], nids)
		# println("after",nids)


		outcome = pick_k_from_array(k-1, nids, nrs, total, cps)


		if outcome == nothing
			continue
		elseif length(outcome) == 0
			continue
		end
	end
	# println("here", total)
	return total

end


function pick_k(k::Int, cps)

	indicies = Array{Int}(undef, 0)
	for i = 1:length(cps)
		push!(indicies, i)
	end

	total = []

	for i = 1:length(indicies)
		result = Array{Int}(undef, 0)
		push!(result, i)
		nids = deepcopy(indicies)

		filter!(x->x≠i, nids)


		pick_k_from_array(k-1, nids, result, total, cps)

	end
	# println(total)
	return total

end

function cm_ids_to_cm(cm_ids, cps)
	cm = []
	# println(cm_ids)
	for i =1 : length(cps)
		if i in cm_ids
			push!(cm, cps[i])
		end
	end

	return myType.cm_id(cm)
end

function find_all_valid_contact_modes(peg::myType.myPeg, socket::myType.mySocket, cs::Array{myType.contact})

	valid_cms = []

	pegs = []

	cps, c_list = find_all_pairs(peg, socket, cs)

	cms = possible_contact_modes(peg, socket, cs)

	# peg_temp = []
	# println(cms)
	c_lists=[]
	# cms = [1, 6]
	for i = 1:length(cms)
		for j =1:length(cms[i])
			cm_ids =cms[i][j]

			# println(cm_ids)

			cm = cm_ids_to_cm(cm_ids, cps)
			# println(cm)
			result, c_list = is_cm_valid(peg, socket, cm)
			# println(result)
			if result
				# println(cm_ids)
				# PyPlot.figure()
				# Draw.draw_sth(peg_temp)
				# Draw.draw_open_sth(socket)

				# println(cm_ids)
				push!(valid_cms, cm_ids)
				# push!(pegs, peg_temp)
				push!(c_lists, c_list)
			end
		end
	end

	return valid_cms, c_lists

end

function find_all_valid_contact_modes(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{myType.contact_id})

	valid_cms = []

	pegs = []

	cps, c_list = find_all_pairs(peg, socket, c_ids)

	cms = possible_contact_modes(peg, socket, c_ids)

	# peg_temp = []
	# println(cms)
	c_lists = []

	# cms = [1, 6]
	for i = 1:length(cms)
		for j =1:length(cms[i])
			cm_ids =cms[i][j]

			# println(cm_ids)

			cm = cm_ids_to_cm(cm_ids, cps)
			# println(cm)
			result, c_list = is_cm_valid(peg, socket, cm)
			# println(result)
			if result
				# println(cm_ids)
				# PyPlot.figure()
				# Draw.draw_sth(peg_temp)
				# Draw.draw_open_sth(socket)

				# println(cm_ids)
				push!(valid_cms, cm_ids)
				# push!(pegs, peg_temp)
				push!(c_lists, c_list)
			end
		end
	end

	return valid_cms, c_lists

end


function find_valid_three_pairs_cms(peg::myType.myPeg, socket::myType.mySocket, c_ids::Array{myType.contact_id})

	valid_cms = []

	pegs = []
	
	cps, c_list = find_all_pairs(peg, socket, c_ids)
	cms = possible_three_pairs_cms(peg, socket, c_ids)

	# peg_temp = []
	# println(cms)
	c_lists = []
	# cms = [1, 6]
	for i = 1:length(cms)
		for j =1:length(cms[i])
			cm_ids =cms[i][j]

			# println(cm_ids)

			cm = cm_ids_to_cm(cm_ids, cps)
			# println(cm)
			result, c_list = is_cm_valid(peg, socket, cm)
			# println(result)
			if result
				# println(cm_ids)
				# PyPlot.figure()
				# Draw.draw_sth(peg_temp)
				# Draw.draw_open_sth(socket)

				# println(cm_ids)
				push!(valid_cms, cm_ids)
				# push!(pegs, peg_temp)
				push!(c_lists, c_list)
			end
		end
	end

	return cps, valid_cms, c_lists
end

function get_title(cps, cms)
	title = []
	n = length(cms)
	for i = 1:n
		title_i = []

		for j = 1: length(cms[i])
			if j == 1
				title_i = "CM: "
			end

			title_i *= "p"*string(cps[cms[i][j]].p_id)*"-"*"e"*string(cps[cms[i][j]].e_id)

			if j == length(cms[i])
				continue
			else
				title_i *= ", "
			end
		end
		push!(title, title_i)
	end

	return  title
end

function cps_to_cp_ids(cps)
	cp_ids = []
	l = length(cps)
	for i=1:l
		push!(cp_ids, [cps[i].p_id, cps[i].p_id])
	end
	return cp_ids

end

function all_combinations(peg, socket, c_ids)
	cps, c_list_useless = find_all_pairs(peg, socket, c_ids)
	cms, c_lists = find_all_valid_contact_modes(peg, socket, c_ids)
	cp_ids = cps_to_cp_ids(cps)
	return cp_ids, cms, c_lists
end




function test_acm(i)
	# peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
    # cs = Construct_Peg_Socket.get_contacts(i)
	# cps, pegs_nothing = find_all_pairs(peg, socket, cs)
	# cms, pegs = find_all_valid_contact_modes(peg, socket, cs)
	#
	# println(cms)
	#
	# title = get_title(cps, cms)
	# # println(title)
	#
	#
	#
	#
	# Draw.draw_all_peg_with_horizontal_line(pegs, socket, title)

	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)
	cps, c_list = find_all_pairs(peg, socket, c_ids)
	cms, c_lists = find_all_valid_contact_modes(peg, socket, c_ids)

	# println(cms)
	# println(c_lists)

	pegs = []
	for i =1:length(c_lists)
		# println(i)
		# println(cms[i])
		# println(c_lists[i][1])
		# for cp in cms[i]
		# 	println(cps[cp])
		# end
		new_peg = Kinematics.get_peg_from_config(peg, c_lists[i][1])
		push!(pegs,new_peg)
	end

	Draw.draw_all_peg(pegs, socket)
	# title = get_title(cps, cms)
	#
	# Draw.draw_all_peg_with_horizontal_line(pegs, socket, title)
end

function test_three_pairs_cms(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Construct_Peg_Socket.get_contact_ids(i)

	cps, cms, c_lists = find_valid_three_pairs_cms(peg, socket, c_ids)

	all_cms = contact_modes_from_three_pairss(cms)
	# println(all_cms)
end

function contact_modes_from_three_pairss(cms::Array)
	n = length(cms)
	all_cms = []
	cms_1 = []
	cms_2 = []
	for cm in cms

		for i = 1:length(cm)-1
			if i == 1
				temp = collect(combinations(cm,i))

				for temp_i in temp
					if !(temp_i in cms_1)
						push!(cms_1, temp_i)
					end
				end
			else
				temp = collect(combinations(cm,i))
				for temp_i in temp
					if !(temp_i in cms_2)
						push!(cms_2, temp_i)
					end
				end
			end
		end
	end
	sort!(cms_1)
	sort!(cms_2)
	for k1 in cms_1
		push!(all_cms, k1)
	end

	for k2 in cms_2
		push!(all_cms, k2)
	end

	for k3 in cms
		push!(all_cms, k3)
	end

	return all_cms
end




function test_acp(i)
	# peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
	# cs = Construct_Peg_Socket.get_contacts(i)
	# cps, pegs = find_all_pairs(peg, socket, cs)
	# println(cps)
	# Draw.draw_all_peg(pegs, socket)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Array{myType.contact_id}(undef, 0)
	push!(c_ids, myType.contact_id(1, 1))
	push!(c_ids, myType.contact_id(2, 2))
	push!(c_ids, myType.contact_id(3, 3))
	push!(c_ids, myType.contact_id(4, 4))
	push!(c_ids, myType.contact_id(5, 5))
	cps, c_list = find_all_pairs(peg, socket, c_ids)
	# println(cps)
	# Draw.draw_all_peg(pegs, socket)
end

function test_cm2(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	c_ids = Array{myType.contact_id}(undef, 0)
	push!(c_ids, myType.contact_id(1, 1))
	push!(c_ids, myType.contact_id(2, 2))
	push!(c_ids, myType.contact_id(3, 3))
	push!(c_ids, myType.contact_id(4, 4))
	push!(c_ids, myType.contact_id(5, 5))
	cps, pegs = find_all_pairs(peg, socket, c_ids)
	# println(cps)
	cm = myType.cm_id([cps[7], cps[8], cps[10]])
	# println(cm)
	result, c_list= is_cm_valid(peg, socket, cm)
	#
	# # dic = distance_dic(peg)
	#
	# println(result)
	#
	# if result
	# 	Draw.draw_all_peg_with_horizontal_line([peg], socket)
	# end
	new_peg = Kinematics.get_peg_from_config(peg, c_list[1])
	# dic = distance_dic(peg)
	Draw.draw_all_peg([new_peg], socket)

	# peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
	# cs = Construct_Peg_Socket.get_contacts(i)
	# cps, c_list = find_all_pairs(peg, socket, cs)
	# println(cps)
	# cm = myType.cm_id([cps[3], cps[6], cps[8]])
	# println(cm)
	# result, c_list= is_cm_valid(peg, socket, cm)
	# new_peg = Kinematics.get_peg_from_config(peg, c_list[1])
	# # dic = distance_dic(peg)
	# Draw.draw_all_peg([new_peg], socket)
	# println(result)

	# if result
	# 	Draw.draw_all_peg_with_horizontal_line([peg], socket)
	# end


end

function test_cm(i)
	peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
	cs = Construct_Peg_Socket.get_contacts(i)
	cps, c_list = find_all_pairs(peg, socket, cs)
	# println(cps)
	cm = myType.cm_id([cps[1]])
	# println(cm)
	result, c_list= is_cm_valid(peg, socket, cm)

	# dic = distance_dic(peg)

	# println(result)
	# if result
	# 	Draw.draw_all_peg_with_horizontal_line([peg], socket)
	# end


end


function test_cp(i)


    # peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
    # cs = Construct_Peg_Socket.get_contacts(i)
	# es = mySocketEdges(socket)
	# cp = contact_pair(1, peg.points[1], es.edges[1])
	#
	# result, new_peg= is_contact_pair_valid(peg, socket, cp)
	# println(result)
	# if result
	# 	Draw.draw_all_peg([new_peg], socket)
	# end
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	cs = Construct_Peg_Socket.get_contacts(i)
	es = mySocketEdges(socket)
	cp = contact_pair(1, peg.points[1], es.edges[1])

	result, config= is_contact_pair_valid(peg, socket, cp)
	peg1 = Kinematics.get_peg_from_config(peg, config[1])
	peg2 = Kinematics.get_peg_from_config(peg, config[length(config)])
	# println(result)
	if result
		Draw.draw_all_peg([peg1,peg2], socket)
	end
end

function test(i)
    # peg, socket = Construct_Peg_Socket.get_peg_and_socket(i)
    # cs = Construct_Peg_Socket.get_contacts(i)
    # println(socket)
    # println(peg)
	peg, socket = Construct_Peg_Socket.read_from_file(i)

	c_ids = Array{myType.contact_id}(undef, 0)
	push!(c_ids, myType.contact_id(1, 1))
	push!(c_ids, myType.contact_id(2, 2))
	push!(c_ids, myType.contact_id(3, 3))
	push!(c_ids, myType.contact_id(4, 4))
	push!(c_ids, myType.contact_id(5, 5))

    es = mySocketEdges(socket)
    center = myPoint2d(peg.points[3].x, peg.points[3].y)
    #result = collision_check(center, peg, socket, 3)
    # cp = contact_pair(1, peg.points[1], es.edges[1])
    # result, new_peg= is_contact_pair_valid(peg, socket, cp)

    cps, c_list = find_all_pairs(peg, socket, c_ids)
	# Constraints.point_on_edge(peg.points[1], es.edges[1])

	# println(cps)
	# Draw.draw_all_peg(pegs, socket)

	# println(result)
	# if result
	# 	Draw.draw_all_peg([new_peg], socket)
	# end



end

function test2()
	find_combination([2,3,1,4], [1,1,1,1], 1, [])
end

function test3(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	# cs = Construct_Peg_Socket.get_contacts(i)
	c_ids = Array{myType.contact_id}(undef, 0)
	push!(c_ids, myType.contact_id(1, 1))
	push!(c_ids, myType.contact_id(2, 2))
	push!(c_ids, myType.contact_id(3, 3))
	push!(c_ids, myType.contact_id(4, 4))
	push!(c_ids, myType.contact_id(5, 5))


	possible_contact_modes(peg, socket, c_ids)


end

end  # module ContactMode
