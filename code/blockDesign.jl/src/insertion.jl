module Insertion


using ..myType
using Compose
using Gadfly
using Colors
using GR

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
import ..Adjust
import ..ContactMode

import PyPlot
import Plots




function compute_signed_angle(v1, v2)

	vn = [0,0,1]


	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return th


end

function test_cp_ag()

	v1 = [-1,1,0]
	v2 = [-1,-1,0]
	th = compute_signed_angle(v1,v2)
	println(th/pi)
end

function all_close(a, b)
	if abs(a - b) < 0.03#change10
		return true
	end
	return false;

end

function pockets_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(-pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

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

	vn = Array{Float64}(undef, length(v))

	for i = 1:length(v)
		vn[i] = -v[i]
	end

	return vn

end


function vector_length(v::Array)
	scale = 0
	for i = 1:length(v)
		scale += v[i]*v[i]
	end

	return sqrt(scale)

end

function inward_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function test_fap(i)

	peg,socket = Construct_Peg_Socket.read_from_file(i)

	hull, pockets = Convexhull.get_convexpocket(socket)
	c = [(1,1,0), (0,1,0), (0,1,1)]
	lw =2.0

	# Visualization.draw_open_sth(hull)
	#
	# for i = 1: length(pockets)
	#
	# 	Visualization.draw_sth(pockets[i],c[2],lw)
	# end
	#
	# Visualization.draw_sth(peg,c[1],lw)


	cps, peg_list=find_all_cp(peg, socket)
	println(cps)

end

function test_cp(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)

	convex = is_convex(socket)
	if !convex
		hull, pockets = Convexhull.get_convexpocket(socket)
	end
	contact_point = myType.Indexed_point(peg.points[3])
	edges = myType.mySocketEdges(socket)
	contact_edge = myType.Indexed_edge(3, socket.points[3], socket.points[4])

	c_p = myType.contact_pair(1, contact_point, contact_edge)
	c_p = myType.contact_pair(1, contact_point, edges.edges[4])

	c = [(1,1,0), (0,1,0), (0,1,1)]
	lw =2.0

	Draw.draw_open_sth(hull)

	for i = 1: length(pockets)

		Draw.draw_sth(pockets[i],c[2],lw)
	end

	Draw.draw_sth(peg,c[1],lw)



	result, c_p, new_peg = is_cp_valid(peg, hull, pockets, socket, c_p)
	# println(result)
	# println(c_p)
	Draw.draw_sth(new_peg,c[3],lw)


end


function inside_socket(peg::myPeg, socket::mySocket,c_p::contact_pair,
	closed = false)
		# use opmization to test if the cp is valid and output valid configu'

		model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

		# unknown variables: the configuration
		JuMP.@variable(model, x)
		JuMP.@variable(model, y)
		JuMP.@variable(model, -pi/6<=th<=pi/6)


		if closed
			n_s = length(socket.edges)
		else
			n_s = length(socket.edges) -1
		end

		for i = 1:n_s
			p_i = socket.points[socket.edges[i].s]
			# get the inward normal
			n = inward_normal(p_i, socket.points[socket.edges[i].e])

			for j = 1:length(peg.points)
				# if j == c_p.c_p.index
				# 	continue
				# end
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

		JuMP.@NLobjective(model, Min, 1 == 1)
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


		for i = 1:n_s
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


function inside_socket(peg::myPeg, socket::mySocket,c_m::contact_mode,
	closed = false)
	# use opmization to test if the MoC is valid and output valid configu'


		  # Draw.draw_sth(peg)
		  # Draw.draw_open_sth(hull)
		  # Draw.draw_open_sth(socket)

		model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

		JuMP.@variable(model, x)
		JuMP.@variable(model, y)
		JuMP.@variable(model, -pi/6<=th<=pi/6)

				if closed
					n_s = length(socket.edges)
				else
					n_s = length(socket.edges) -1
				end


			for i = 1:n_s
				p_i = socket.points[socket.edges[i].s]
				# get the inward normal
				n = inward_normal(p_i, socket.points[socket.edges[i].e])

				for j = 1:length(peg.points)
					# if point_in_contact_mode(j, c_m)
					# 	continue
					# end
					# for each potential contact point on peg
					# add constraint of the positive
					cp_j = peg.points[j]

					JuMP.@NLconstraint(model,
						n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
						n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0.001)

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

			JuMP.@NLobjective(model, Min, 1 == 1)
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
			for i = 1:n_s
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





end


function is_cp_valid(peg::myPeg, socket::mySocket, c_p::contact_pair,
	closed=false)

	# test the cp for convex joint

	#PyPlot.figure()


	result=false

	count = 0

	all_vio = Array{Array{Int}}(undef, 0)



	vx,vy,vth = inside_socket(peg, socket, c_p, closed)


	if vx==nothing
		#PyPlot.close()


		return result, c_p, peg
	end

	config = configSE2(vx,vy,vth)
 	new_peg = Kinematics.get_peg_from_config(peg,config)
 	# Draw.draw_mlti_sth(socket_list, new_peg)
	# Draw.draw_sth(peg)


	result = !result
	#
	# if result == true
	#
	#
	# 	Draw.draw_sth(new_peg,(0,1,1),2.0)
	#
	# else
	# 	println("the cp is not valid")
	# end

	#PyPlot.pause(1)

	#PyPlot.close()
	return result, c_p, new_peg



end

function is_cp_valid(peg::myPeg, hull::mySocket, pockets::Array{myType.mySocket},
	socket::mySocket, c_p::contact_pair, closed = false)
	#test cp for concave joint
	#for the concave joint, the peg should in the hull but out the pockets



	#PyPlot.figure()

	socket_list = Array{mySocket}(undef, 0)

	push!(socket_list, hull)
	for i=1:length(pockets)
		push!(socket_list,pockets[i])
	end

	result=false

	count = 0

	all_vio = Array{Array{Int}}(undef, 0)



	vx,vy,vth = inside_convexhull(peg, hull, socket, c_p, closed)


	if vx==nothing
		#PyPlot.close()


		return result, c_p, peg
	end

	config = configSE2(vx,vy,vth)
 	new_peg = Kinematics.get_peg_from_config(peg,config)
 	# Draw.draw_mlti_sth(socket_list, new_peg)
	# Draw.draw_sth(peg)


	vio,vio_id = test_vio(new_peg, pockets, socket, c_p)

	if vio == true

		for i = length(vio_id)
			push!(all_vio, vio_id[i])
		end
	end





	while vio == true
		vx,vy,vth = outside_pockets(new_peg, hull, socket, c_p, all_vio, closed)




		# if vx==nothing
		# 	#PyPlot.close()
		# 	println("here")
		# 	return result, c_p, new_peg
		# end


		config = configSE2(vx,vy,vth)

		new_peg = Kinematics.get_peg_from_config(new_peg,config)

	 	vio,vio_id = test_vio(new_peg, pockets, socket, c_p)

		if vio == true

			for i = length(vio_id)
				push!(all_vio, vio_id[i])
			end
		end

		count += 1

		if count >=30
			#PyPlot.close()

			return result, c_p, peg
		end


	end

	result = !result
	#
	# if result == true
	#
	#
	# 	Draw.draw_sth(new_peg,(0,1,1),2.0)
	#
	# else
	# 	println("the cp is not valid")
	# end

	#PyPlot.pause(1)

	#PyPlot.close()
	return result, c_p, new_peg



end


function is_mode_valid(peg::myPeg, socket::mySocket, c_m::contact_mode, closed = false)
	# test the MoC for convex joint


	# socket_list = Array{mySocket}(undef, 0)
	#
	# push!(socket_list, hull)
	# for i=1:length(pockets)
	# 	push!(socket_list,pockets[i])
	# end
	#
	# PyPlot.figure()



	#cf_list = Array{myType.configSE2}(undef, 0)

	result=false

	count = 0

	all_vio = Array{Array{Int}}(undef, 0)



	vx,vy,vth = inside_socket(peg, socket, c_m, closed)


	if vx==nothing


		return result, c_m, peg
	end

	config = configSE2(vx,vy,vth)
	new_peg = Kinematics.get_peg_from_config(peg,config)


	# Draw.draw_open_sth(socket)
	# Draw.draw_sth(new_peg)


	result = !result

	# if result == true
	#
	# 	Draw.draw_sth(new_peg,(0,1,1),2.0)
	#
	# else
	# 	println("the cp is not valid")
	# end

	return result, c_m, new_peg



end


function is_mode_valid(peg::myPeg, hull::mySocket, pockets::Array{myType.mySocket},
	socket::mySocket, c_m::contact_mode, closed = false)
	#test MoC for concave joint
	#for the concave joint, the peg should in the hull but out the pockets

	# socket_list = Array{mySocket}(undef, 0)
	#
	# push!(socket_list, hull)
	# for i=1:length(pockets)
	# 	push!(socket_list,pockets[i])
	# end
	#
	# PyPlot.figure()



	#cf_list = Array{myType.configSE2}(undef, 0)

	result=false

	count = 0

	all_vio = Array{Array{Int}}(undef, 0)



	vx,vy,vth = inside_convexhull(peg, hull, socket, c_m, closed)


	if vx==nothing


		return result, c_m, peg
	end

	config = configSE2(vx,vy,vth)
	new_peg = Kinematics.get_peg_from_config(peg,config)


	# Draw.draw_open_sth(socket)
	# Draw.draw_sth(new_peg)

	vio,vio_id = test_vio(new_peg, pockets, socket, c_m)

	if vio == true

		for i = length(vio_id)
			push!(all_vio, vio_id[i])
		end
	end





	while vio == true && count <=10
		vx,vy,vth = outside_pockets(new_peg, hull, socket, c_m, all_vio, closed)




		# if vx==nothing
		#
		# 	return result, c_m, new_peg
		# end


		config = configSE2(vx,vy,vth)

		new_peg = Kinematics.get_peg_from_config(new_peg,config)

	 	vio,vio_id = test_vio(new_peg, pockets, socket, c_m)
		#println(vio_id)


		if vio == true

			for i = length(vio_id)
				push!(all_vio, vio_id[i])
			end
		end

		count += 1

		# if count >=3
		#
		# 	return result, c_m, peg
		# end


	end


	if vio ==  false
		result = !result
	end

	# if result == true
	#
	# 	Draw.draw_sth(new_peg,(0,1,1),2.0)
	#
	# else
	# 	println("the cp is not valid")
	# end

	return result, c_m, new_peg



end




function inside_convexhull(peg::myPeg, hull::mySocket, socket::mySocket,
	c_p::contact_pair, closed = false)

	# It's similar with inside socket


	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, -pi/6<=th<=pi/6)

	if closed
		n_h = length(hull.edges)
	else
		n_h = length(hull.edges) -1
	end

	for i = 1:n_h
		p_i = socket.points[hull.edges[i].s]
		n = inward_normal(p_i, socket.points[hull.edges[i].e])

		for j = 1:length(peg.points)
			if j == c_p.c_p.index && i ==c_p.c_e.index
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			JuMP.@NLconstraint(model,
				n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
				n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0);#change1


		end
	end




	p_i = c_p.c_e.s
	n = inward_normal(p_i, c_p.c_e.e)

	p = peg.points[c_p.c_p.index]

	min_x = min(p_i.x, c_p.c_e.e.x)
	max_x = max(p_i.x, c_p.c_e.e.x)
	min_y = min(p_i.y, c_p.c_e.e.y)
	max_y = max(p_i.y, c_p.c_e.e.y)

	#println(p,",",min_x,",",max_x,",",min_y,",",max_y)

	# last add the equality constraint of the contact pair
	JuMP.@NLconstraint(model,
		n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
		n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) == 0);

	JuMP.@NLconstraint(model,
		min_x<=(cos(th) * p.x - sin(th) * p.y + x)<=max_x
		);

	JuMP.@NLconstraint(model,
		min_y<=(sin(th) * p.y + cos(th) * p.y + y)<=max_y
		);






	# minimize the absolute value of theta

	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	vio = false
	#println(c_p)

	config = configSE2(vx,vy,vth)
	new_peg = Kinematics.get_peg_from_config(peg,config)
	#Draw.draw_sth(new_peg)






	for i = 1:n_h
		p_i = socket.points[hull.edges[i].s]

		# get the inward normal
		n = inward_normal(p_i, socket.points[hull.edges[i].e])

		for j = 1:length(peg.points)
			if j == c_p.c_p.index
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]


			if n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
				n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) < -0.001
				vio = true
			end

		end
	end

	p_i = c_p.c_e.s
	n = inward_normal(p_i, c_p.c_e.e)

	p = peg.points[c_p.c_p.index]
	# last add the equality constraint of the contact pair
	# println("eq: ", n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
	# 	n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y) );

	#println(n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
	#		n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y))

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

function inside_convexhull(peg::myPeg, hull::mySocket, socket::mySocket,
	c_m::contact_mode, closed = false)
	  # Draw.draw_sth(peg)
	  # Draw.draw_open_sth(hull)
	  # Draw.draw_open_sth(socket)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, -pi/6<=th<=pi/6)


	if closed
		n_h = length(hull.edges)
	else
		n_h = length(hull.edges) -1
	end

	# NEED TO MODIFy



	for i = 1:n_h
		p_i = socket.points[hull.edges[i].s]
		n = inward_normal(p_i, socket.points[hull.edges[i].e])

		for j = 1:length(peg.points)
			if point_in_contact_mode(j, c_m)
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			JuMP.@NLconstraint(model,
				n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
				n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0.001);#change2


		end
	end




	for i = 1:length(c_m.cps)
		#println(c_m.cps[i])

		p_i = c_m.cps[i].c_e.s
		n = inward_normal(p_i, c_m.cps[i].c_e.e)

		# n2 = [c_m.cps[i].c_e.e.x-p_i.x, c_m.cps[i].c_e.e.y-p_i.y]
		#
		# n2 = LinearAlgebra.normalize(n2)
		#
		#
		# d = sqrt(n[1]^2+n[2]^2)


		p = peg.points[c_m.cps[i].c_p.index]

		min_x = min(p_i.x, c_m.cps[i].c_e.e.x)
		max_x = max(p_i.x, c_m.cps[i].c_e.e.x)
		min_y = min(p_i.y, c_m.cps[i].c_e.e.y)
		max_y = max(p_i.y, c_m.cps[i].c_e.e.y)

		#println(p,",",min_x,",",max_x,",",min_y,",",max_y)

		# last add the equality constraint of the contact pair
		JuMP.@NLconstraint(model,
			n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
			n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) >= 0);

		JuMP.@NLconstraint(model,
			n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
			n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) <= 0);

		# JuMP.@NLconstraint(model,
		# 	n2[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
		# 	n2[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) >= 0);
		#
		# JuMP.@NLconstraint(model,
		# 	n2[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
		# 	n2[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) <= d);
		#
		JuMP.@NLconstraint(model,
			min_x<=(cos(th) * p.x - sin(th) * p.y + x)<=max_x
			);
		# #
		# JuMP.@NLconstraint(model,
		# 	min_y<=(sin(th) * p.y + cos(th) * p.y + y)<=max_y
		# 	);

	end





		# minimize the absolute value of theta

	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	vio = false
	#println(c_p)
	config = configSE2(vx,vy,vth)
	new_peg = Kinematics.get_peg_from_config(peg,config)

	# Visualization.draw_sth(new_peg)





	for i = 1:n_h
		p_i = socket.points[hull.edges[i].s]
		# get the inward normal
		n = inward_normal(p_i, socket.points[hull.edges[i].e])

		for j = 1:length(peg.points)
			if point_in_contact_mode(j, c_m)
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			# println(n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
			# 	n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) )
			if n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
				n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) < 0#change3
				vio = true
			end

		end
	end

	for i = 1:length(c_m.cps)
		p_i = c_m.cps[i].c_e.s
		n = inward_normal(p_i, c_m.cps[i].c_e.e)


		p = peg.points[c_m.cps[i].c_p.index]
	# last add the equality constraint of the contact pair
	# println("eq: ", n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
	# 	n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y) );

	#println(n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
	#		n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y))

		if !all_close(n[1] * ((cos(vth) * p.x - sin(vth) * p.y + vx) - p_i.x) +
				n[2] * ((sin(vth) * p.x + cos(vth) * p.y + vy) - p_i.y), 0)
				vio = true
		end
	end



	if status == :Optimal && !vio
		return vx, vy, vth
	end

	return nothing, nothing, nothing





end


function test_vio(peg::myType.myPeg, pockets::Array{myType.mySocket},socket::myType.mySocket,
	c_p::myType.contact_pair)

	# test if the point of the peg in the socket, if yes, output its index



	vio = false

	vio_id = Array{Array{Int}}(undef, 0)

	for i=1:length(pockets)

		for j = 1:length(peg.points)
			if j == c_p.c_p.index
				continue
			end


			p_j =peg.points[j]

			dis_set = []

			for l = 1:length(pockets[i].edges)-1

				p_l = socket.points[pockets[i].edges[l].s]

				n = pockets_normal(p_l, socket.points[pockets[i].edges[l].e])



				v = [p_j.x-p_l.x, p_j.y-p_l.y]
				dis = n[1]*v[1]+n[2]*v[2]
				#if dis <=0
				#	break
				#end

				push!(dis_set, [dis, j, pockets[i].edges[l].index])
			end

			sort!(dis_set, by = x -> x[1])

			#println(dis_set)


			if dis_set[1][1]<= 0.01#change4-0.01
				continue
			end
			#println(dis_set)

			vio = true


			if !([dis_set[1][2],dis_set[1][3]] in vio_id)
				push!(vio_id, [dis_set[1][2],dis_set[1][3]])
			end

		end
	end

	return vio, vio_id


end

function test_vio(peg::myType.myPeg, pockets::Array{myType.mySocket},
	socket::myType.mySocket, c_m::myType.contact_mode)

	# same with test_vio for cp

	vio = false

	vio_id = Array{Array{Int}}(undef, 0)

	for i=1:length(pockets)

		for j = 1:length(peg.points)
			if point_in_contact_mode(j, c_m)
				continue
			end


			p_j =peg.points[j]

			dis_set = []

			for l = 1:length(pockets[i].edges)-1

				p_l = socket.points[pockets[i].edges[l].s]

				n = pockets_normal(p_l, socket.points[pockets[i].edges[l].e])



				v = [p_j.x-p_l.x, p_j.y-p_l.y]
				dis = n[1]*v[1]+n[2]*v[2]
				#if dis <=0
				#	break
				#end

				push!(dis_set, [dis, j, pockets[i].edges[l].index])
			end

			sort!(dis_set, by = x -> x[1])
			#println(dis_set)

			if dis_set[1][1]<= 0.03
				continue
			end

			vio = true
			# println(vio)
			 #println(dis_set)

			if !([dis_set[1][2],dis_set[1][3]] in vio_id)
				push!(vio_id, [dis_set[1][2],dis_set[1][3]])
			end

		end
	end

	# println("here")
	# println(vio)
	return vio, vio_id


end



function outside_pockets(peg::myType.myPeg, hull::myType.mySocket,
	socket::myType.mySocket, c_p::myType.contact_pair,
	vio_id, closed=false)

	# use optimization to make peg out the pokects

	# Draw.draw_open_sth(socket)
	# Draw.draw_sth(peg)
	#
	# println(vio_id)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, -pi/6<=th<=pi/6)

	if closed
		n_h = length(hull.edges)
	else
		n_h = length(hull.edges) -1
	end

	for i = 1:n_h

		p_i = socket.points[hull.edges[i].s]
		n = inward_normal(p_i, socket.points[hull.edges[i].e])

		for j = 1:length(peg.points)


			if j == c_p.c_p.index #|| j == vio_id[1]
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			JuMP.@NLconstraint(model,
				n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
				n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0);#change7


		end
	end



	for i =1:length(vio_id)


		ps_vio = socket.points[socket.edges[vio_id[i][2]].s]
		nv = inward_normal(ps_vio, socket.points[socket.edges[vio_id[i][2]].e])
		p_vio = peg.points[vio_id[i][1]]

		JuMP.@NLconstraint(model,
			nv[1] * ((cos(th) * p_vio.x - sin(th) * p_vio.y + x) - ps_vio.x) +
			nv[2] * ((sin(th) * p_vio.x + cos(th) * p_vio.y + y) - ps_vio.y) >= 0);#change6
	end





	p_i = c_p.c_e.s
	n = inward_normal(p_i, c_p.c_e.e)

	p = peg.points[c_p.c_p.index]
	# last add the equality constraint of the contact pair

	min_x = min(p_i.x, c_p.c_e.e.x)
	max_x = max(p_i.x, c_p.c_e.e.x)
	min_y = min(p_i.y, c_p.c_e.e.y)
	max_y = max(p_i.y, c_p.c_e.e.y)

			# last add the equality constraint of the contact pair
	JuMP.@NLconstraint(model,
		n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
		n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) == 0);

	JuMP.@NLconstraint(model,
		min_x<=(cos(th) * p.x - sin(th) * p.y + x)<=max_x
		);

	JuMP.@NLconstraint(model,
		min_y<=(sin(th) * p.y + cos(th) * p.y + y)<=max_y
		);


		# minimize the absolute value of theta

		JuMP.@NLobjective(model, Min, 1==1)


		status = JuMP.solve(model)


		vx = JuMP.getvalue(x)
		vy = JuMP.getvalue(y)
		vth = JuMP.getvalue(th)


		config = configSE2(vx,vy,vth)

		new_peg = Kinematics.get_peg_from_config(peg,config)

		# Draw.draw_sth(new_peg)
		# Draw.draw_open_sth(socket)


		return vx, vy, vth



end



function outside_pockets(peg::myType.myPeg, hull::myType.mySocket,
	socket::myType.mySocket, c_m::myType.contact_mode,
	vio_id, closed = false)

	#It's the same with outside_pockets function for cp

	#Draw.draw_sth(peg)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, -pi/6<=th<=pi/6)

	if closed
		n_h = length(hull.edges)
	else
		n_h = length(hull.edges) -1
	end

	for i = 1:n_h

		p_i = socket.points[hull.edges[i].s]
		n = inward_normal(p_i, socket.points[hull.edges[i].e])

		for j = 1:length(peg.points)


			if point_in_contact_mode(j, c_m) #|| j == vio_id[1]
				continue
			end
			# for each potential contact point on peg
			# add constraint of the positive
			cp_j = peg.points[j]

			JuMP.@NLconstraint(model,
				n[1] * ((cos(th) * cp_j.x - sin(th) * cp_j.y + x) - p_i.x) +
				n[2] * ((sin(th) * cp_j.x + cos(th) * cp_j.y + y) - p_i.y) >= 0.01);


		end
	end



	for i =1:length(vio_id)


		ps_vio = socket.points[socket.edges[vio_id[i][2]].s]
		nv = inward_normal(ps_vio, socket.points[socket.edges[vio_id[i][2]].e])
		p_vio = peg.points[vio_id[i][1]]

		JuMP.@NLconstraint(model,
			nv[1] * ((cos(th) * p_vio.x - sin(th) * p_vio.y + x) - ps_vio.x) +
			nv[2] * ((sin(th) * p_vio.x + cos(th) * p_vio.y + y) - ps_vio.y) >= 0.01);
	end



	for i = 1:length(c_m.cps)
		p_i = c_m.cps[i].c_e.s
		n = inward_normal(p_i, c_m.cps[i].c_e.e)


		p = peg.points[c_m.cps[i].c_p.index]

		min_x = min(p_i.x, c_m.cps[i].c_e.e.x)
		max_x = max(p_i.x, c_m.cps[i].c_e.e.x)
		min_y = min(p_i.y, c_m.cps[i].c_e.e.y)
		max_y = max(p_i.y, c_m.cps[i].c_e.e.y)

		#println(p,",",min_x,",",max_x,",",min_y,",",max_y)

		# last add the equality constraint of the contact pair
		JuMP.@NLconstraint(model,
			n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
			n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) >= 0);

		JuMP.@NLconstraint(model,
			n[1] * ((cos(th) * p.x - sin(th) * p.y + x) - p_i.x) +
			n[2] * ((sin(th) * p.x + cos(th) * p.y + y) - p_i.y) <= 0);

		JuMP.@NLconstraint(model,
			min_x<=(cos(th) * p.x - sin(th) * p.y + x)<=max_x
			);

		JuMP.@NLconstraint(model,
			min_y<=(sin(th) * p.y + cos(th) * p.y + y)<=max_y
			);

	end


		# minimize the absolute value of theta

		JuMP.@NLobjective(model, Min, 1==1)


		status = JuMP.solve(model)


		vx = JuMP.getvalue(x)
		vy = JuMP.getvalue(y)
		vth = JuMP.getvalue(th)


		config = configSE2(vx,vy,vth)

		new_peg = Kinematics.get_peg_from_config(peg,config)

		# Draw.draw_sth(new_peg)
		# Draw.draw_open_sth(socket)


	return vx, vy, vth


end



function find_all_cp(peg::myType.myPeg, socket::myType.mySocket, closed=false, convex=false)

	# find all valid cp for the peg and socket

	peg_list = Array{myType.myPeg}(undef, 0)


	if !convex
		hull, pockets = Convexhull.get_convexpocket(socket)
	end

	lp = length(peg.points)
	if closed
		le = length(socket.edges)
	else
		le = length(socket.edges)- 1
	end

	se = myType.mySocketEdges(socket)

	ps = []

	for i = 1:lp
		for j = 1:le
			# try each possible combination
			# test pair

			c_p = myType.contact_pair((i-1)*lp + j, peg.points[i], se.edges[j])
			if !convex
				result, c_p, new_peg = is_cp_valid(peg, hull, pockets, socket, c_p, closed)
			else
				result, c_p, new_peg = is_cp_valid(peg, socket, c_p, closed)
			end


			# println(result)
			# println([i,j], c_p)
			if result == false
				continue
			end


			# pts = Kinematics.compute_points_from_config(peg, myType.configSE2(vx, vy, vth))


			# Visualization.display_peg_socket(socket, pts)
			# sleep(1)

			# otherwise, it exists,
			push!(ps, [i, j])
			push!(peg_list, new_peg)
		end
	end

	return ps, peg_list

	# println(ps)


end


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





function validate_mode(list, peg::myType.myPeg, socket::myType.mySocket,
	c_p_s, closed=false, convex = false)

	# out put all valid pegs with the MoC in list


	peg_list=Array{myType.myPeg}(undef, 0)

	se = myType.mySocketEdges(socket)
	# println("enter: ", length(list))
	# list is a mode
	# not a list of modes
	if !convex
		hull, pockets = Convexhull.get_convexpocket(socket)
	end


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

		if !convex
			result, c_m, new_peg = is_mode_valid(peg, hull, pockets, socket, c_m, closed)
		else
			result, c_m, new_peg = is_mode_valid(peg, socket, c_m, closed)
		end
		#println(mode)


		#println(result)


		if result == false
			# remove this from the list
			push!(invalid, mode)

		else
			push!(peg_list, new_peg)

		end





	end

	for mode in invalid
		filter!(x->x≠mode, list)
	end

	return peg_list

	# println("exiting: ", length(list))

end



function all_combinations(peg::myType.myPeg, socket::myType.mySocket, closed=false, convex=false)
	# for each contact point, form with each of the edge
	# also, all possible number of pairs
	# how to loop?

	# first loop over the possible number of pairs

	# contact pairs

	model_list = Array{Array}(undef, 0)

	all_peg = Array{myType.myPeg}(undef, 0)


	c_p_s, pl = find_all_cp(peg, socket, closed, convex)




	 for i = 1:length(pl)
		push!(all_peg, pl[i])
	 end
	# println(length(c_p_s))
	#println(c_p_s)
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
		#println(total)

		peg_list = validate_mode(total, peg, socket, c_p_s, closed, convex)

		if !isempty(total)

			push!(model_list, total)


			for j = 1:length(peg_list)
				push!(all_peg, peg_list[j])
			end
		end


		#println(total)


	end

	return all_peg, c_p_s, model_list




	# for every result in the mode
	# test if the mode is possible

	# total = pick_k(5, c_p_s)

end

function test_fam(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
	all_peg, c_p_s, model_list = all_combinations(peg, socket)

	Draw.draw_all_peg(all_peg, peg, socket)






end

function test_mode(i)

	peg,socket = Construct_Peg_Socket.read_from_file(i)


	#peg = Adjust.get_adj_peg(peg, socket, 0.08)
	convex = is_convex(socket)
	if !convex
		hull, pockets = Convexhull.get_convexpocket(socket)
	end

	closed = false
	se = myType.mySocketEdges(socket)

	c_p_s, peg_list = find_all_cp(peg, socket, closed, convex)
	# println(c_p_s)

	contact_pairs = Array{myType.contact_pair}(undef, 0)




	push!(contact_pairs, myType.contact_pair(1, peg.points[c_p_s[1][1]], se.edges[c_p_s[1][2]]))
	push!(contact_pairs, myType.contact_pair(10, peg.points[c_p_s[10][1]], se.edges[c_p_s[10][2]]))

	c_m = myType.contact_mode(contact_pairs)
	if !convex
		result, c_m, new_peg = is_mode_valid(peg, hull, pockets, socket, c_m, closed)
	else
		result, c_m, new_peg = is_mode_valid(peg, socket, c_m, closed)
	end
	config = myType.configSE2(-0.6, 0.6, 0.1)
	new_peg = Kinematics.get_peg_from_config(new_peg, config)

	Draw.draw_sth(new_peg, (0,0,1), 3.)
	Draw.draw_open_sth(socket, (0,0,0), 3.)
	# println("result:", result)

end


function move_peg(peg::myType.myPeg)
	n = length(peg.points)
	for i = 1:n
		peg.points[i].x-=0.6
		peg.points[i].y+=0.6
	end

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


	# num_diff = 0
	# for i = 1:l1
	# 	current = cm1.cps[i]
	# 	for j = 1:l2
	# 		if cm2.cps[j].c_p.index < current.c_p.index
	# 			continue
	# 		elseif cm2.cps[j].c_p.index == current.c_p.index &&
	# 			cm2.cps[j].c_e.index < current.c_e.index
	# 			continue
	# 		elseif cm2.cps[j].c_p.index == current.c_p.index &&
	# 			cm2.cps[j].c_e.index == current.c_e.index
	# 			break
	# 		elseif cm2.cps[j].c_p.index == current.c_p.index &&
	# 			cm2.cps[j].c_e.index > current.c_e.index
	# 			num_diff += 1
	# 			break
	# 		elseif cm2.cps[j].c_p.index > current.c_p.index
	# 			num_diff += 1
	# 			break
	# 		end
	# 	end
	# end
	#
	# if num_diff > 1
	# 	return false
	# end

	return true

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

	if scale == 0
		return v
	end

	factor = target_scale / scale
	#println(factor)

	u = deepcopy(v)

	for i = 1:length(v)

		u[i] = u[i] * abs(factor)

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
		# println("cm1 and cm2 are not neighbors")
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

	#println(px,py)

	c_p = myType.Indexed_point(contact_point.index, px, py)


	dist = dist_cp_ce(c_p, contact_edge, peg, socket)

	#println(dist)

	if length(cm1.cps) < length(cm2.cps)

		if all_close(dist, 0)

			# println(c_p)
			# println(contact_edge)
			p_id = Array{Int}(undef, 0)
			for i = 1:length(cm1.cps)
				push!(p_id, cm1.cps[i].c_p.index)
			end
			if c_p.index in p_id
				return true
			else
				return trans_when_dis0(peg, socket, cm1, cm2, c_p, contact_edge, config, F)
			end

		end


		is_cp = false

		F_all = [0,0]

		cp_id = 0

		for i = 1:length(cm1.cps)

			if contact_point.index == cm1.cps[i].c_p.index
				is_cp = true
				cp_id = i
				continue
			end


			F_i = deepcopy(F)

			# for j = 1:length(cm1.cps)
			# 	if j == i || contact_point.index == cm1.cps[j].c_p.index
			# 		continue
			# 	end
			# 	F_pj, Fnj = force_between_pair(cm1.cps[j], peg, socket, config, F)
			# 	if Fnj >0
			# 		F_j += F_pj
			# 	end
			# end



			F_pi, Fni= force_between_pair(cm1.cps[i], peg, F_i)

			#println(vector_length(F_pi))
			if Fni>0
				F_all += F_pi
			end
		end

		if !is_cp
			F_p = [F_all[1]+F[1], F_all[2]+F[2]]
		else

			F_all = [F_all[1]+F[1], F_all[2]+F[2], F[3]]

			F_p, Fn = force_between_pair(cm1.cps[cp_id], peg, socket, F_all)


		end

		#println(vector_length([F[1], F[2]]),",", F_p)
		# if all_close(dist,0)
		# 	scale = 1
		# else
		# 	scale =dist/2
		# end
		# println(F_p)
		# println(Fn)
		F_s = rescale_vector(F_p, dist/2)


		ncp = myType.Indexed_point(c_p.index, c_p.x + F_s[1],
								c_p.y + F_s[2])


		n_dist = dist_cp_ce(ncp, contact_edge, peg, socket)

		o_dist = 0

		if is_cp
			o_dist = dist_cp_ce(ncp, cm1.cps[cp_id].c_e, peg, socket)

		end


		if n_dist >=dist || !all_close(o_dist,0)
			return false
		end


		# if (n_dist > dist)||all_close(n_dist,dist)
		# 	# increased distance
		# 	return false
		# end




	elseif length(cm1.cps) > length(cm2.cps)
		F_all = [0,0]

		cp_id = 0
		# in_two_cp = false
		#
		# for i = 1:length(cm2.cps)
		#
		# 	if contact_point.index == cm2.cps[i].c_p.index
		# 		in_two_cp = true
		# 	end
		# end

		#println(in_two_cp)
		for i = 1:length(cm1.cps)

			# if in_two_cp == true && contact_edge.index == cm1.cps[i].c_e.index
			# 	continue
			# end

			if contact_point.index == cm1.cps[i].c_p.index
				cp_id = i
				continue
			end

			F_j = [0,0]
			F_i = deepcopy(F)

			# for j = 1:length(cm1.cps)
			# 	if j == i || cm1.cps[j].c_p.index == contact_point.index
			# 		continue
			# 	end
			# 	F_pj, Fnj = force_between_pair(cm1.cps[j], peg, socket, config, F)
			#
			# 	if Fnj>0
			# 		F_j += F_pj
			# 	end
			#
			# end

			F_i = [F[1]+F_j[1], F[2]+F_j[2], F[3]]

			F_pi, Fni = force_between_pair(cm1.cps[i], peg, socket, F_i)
			#println(vector_length(F_pi))
			if Fni>0
				F_all += F_pi
			end

		end

		F_all = [F_all[1]+F[1], F_all[2]+F[2], F[3]]
		#
		F_p, Fn = force_between_pair(cm1.cps[cp_id], peg, socket, F_all)

		if Fn > 0
			return false
		end


		# if all_close(dist,0)
		# 	scale = 1
		# else
		# 	scale =dist/2
		# end

		# F_s = rescale_vector(F_p, 100000)
		#
		# ncp = myType.Indexed_point(c_p.index, c_p.x + F_s[1],
		# 						c_p.y + F_s[2])
		#
		#
		# n_dist = dist_cp_ce(ncp, contact_edge, peg, socket)
		#
		#
		#
		# if all_close(n_dist,dist)
		#
		# 	return false
		# end

	elseif length(cm1.cps) == length(cm2.cps)
		# println(cm1.cps)
		# println(cm2.cps)
		# println("no such transfer")
		return false
	end




	#
	#
	# for i = 1:length(cm1.cps)
	#
	# 	F_p = force_between_pair(cm1.cps[i], peg, socket, config, F)
	# 	println(c_p)
	#
	# 	println(F_p)
	#
	# 	# this is the force at contact point
	# 	#
	#
	# 	# need to normalize F_p so that it does not overshoot
	#
	# 	# this F_s is the scaled force
	# 	F_s = rescale_vector(F_p, dist / 2)
	#
	# 	println(F_s)
	# 	ncp = myType.Indexed_point(c_p.index, c_p.x + F_s[1],
	# 							c_p.y + F_s[2])
	#
	# 	n_dist = dist_cp_ce(ncp, contact_edge, peg, socket)
	# 	println(ncp)
	# 	println(n_dist)
	# 	println(dist)
	# 	# if !is_valid_force(F_p, cm1.cps[i], c_p, contact_edge)
	# 	# 	return false
	# 	# end
	#
	# 	if (n_dist > dist)||all_close(n_dist,dist)
	# 		# increased distance
	# 		return false
	# 	end
	#
	# end


	return true


end

function trans_when_dis0_pp(peg::myType.myPeg, socket::myType.mySocket, cm1::myType.contact_mode, cm2::myType.contact_mode,
							c_p::myType.Indexed_point, c_e::myType.Indexed_edge, config::myType.configSE2, F::Array)

	# transition between two MoC with point - point situation


	id = c_p.index
	px = c_p.x
	py = c_p.y
	new_peg = Kinematics.get_peg_from_config(peg, config)
	p_c = compute_centroid(new_peg)

	v = [px-p_c[1], py-p_c[2]]

	F_all = [0,0]

	id_set = Array{Int}(undef, 0)

	push!(id_set, c_e.index)


	for i = 1:length(cm1.cps)

		if id == cm1.cps[i].c_p.index
			push!(id_set, cm1.cps[i].c_e.index)
			continue
		end

		F_j = [0,0]
		F_i = deepcopy(F)

		# for j = 1:length(cm1.cps)
		# 	if j == i||id == cm1.cps[j].c_p.index
		# 		continue
		# 	end
		#
		# 	F_pj, Fnj = force_between_pair(cm1.cps[j], peg, socket, config, F)
		#
		# 	if Fnj > 0
		# 		F_j += F_pj
		# 	end
		#
		# end

		F_i = [F[1]+F_j[1], F[2]+F_j[2], F[3]]

		F_pi, Fni = force_between_pair(cm1.cps[i], peg, socket, F_i)

		if Fni > 0
			F_all += F_pi
		end

	end

	F_all = [F_all[1]+F[1], F_all[2]+F[2]]

	e_id =c_e.index

	sort!(id_set)

	vn = Array{Array{Float64}}(undef, 0)
	for i in id_set
		next = (i+1)%length(socket.points)
		if next == 0
			next = length(socket.points)
		end
		x = socket.points[i].x - socket.points[next].x
		y = socket.points[i].y - socket.points[next].y
		push!(vn, [x, y, 0])
	end

	vn[2]=-vn[2]


	th = compute_signed_angle([F_all[1], F_all[2], 0], vn[2])


	th1 = compute_signed_angle(vn[1], vn[2])/2

	if th < 0
		th += 2pi
	end

	if th1 < 0
		th1 += 2pi
	end


	if th >= th1 && e_id == id_set[1]
		return true
	elseif th <= th1 && e_id == id_set[2]
		return true
	end

	return false

end

function trans_when_dis0(peg::myType.myPeg, socket::myType.mySocket, cm1::myType.contact_mode, cm2::myType.contact_mode,
							c_p::myType.Indexed_point, c_e::myType.Indexed_edge, config::myType.configSE2, F::Array)

	# transition between two MoC with the distance = 0
	id = c_p.index
	px = c_p.x
	py = c_p.y
	new_peg = Kinematics.get_peg_from_config(peg, config)

	F_all = [0,0]


	for i = 1:length(cm1.cps)


		F_j = [0,0]
		F_i = deepcopy(F)

		# for j = 1:length(cm1.cps)
		# 	if j == i||id == cm1.cps[j].c_p.index
		# 		continue
		# 	end
		#
		# 	F_pj, Fnj = force_between_pair(cm1.cps[j], peg, socket, config, F)
		#
		# 	if Fnj > 0
		# 		F_j += F_pj
		# 	end
		#
		# end

		F_i = [F[1]+F_j[1], F[2]+F_j[2], F[3]]

		F_pi, Fni = force_between_pair(cm1.cps[i], peg, socket, F_i)

		if Fni > 0
			F_all += F_pi
		end

	end

	F_all = [F_all[1]+F[1], F_all[2]+F[2]]


	# next = (id+1)%length(socket.points)
	# if next == 0
	# 	next = length(socket.points)
	# end

	# vn = [socket.points[next].x - socket.points[id].x,
	# 	socket.points[next].y - socket.points[id].y,
	# 	0]
	vn = [c_e.e.x - c_e.s.x, c_e.e.y - c_e.s.y, 0]

	th = compute_signed_angle([F_all[1], F_all[2], 0], vn)
	# println(c_p.index, ",", c_e.index)
	# println(F_all,",", vn)
 	# println(th)

	if th < 0
		return false
	else
		return true
	end

	return false

end

function p_p_direction(cm1::myType.contact_mode, cm2::myType.contact_mode,
						peg::myType.myPeg, socket::myType.mySocket,
						config::myType.configSE2, F::Array, id::Int)
	# point-point situation

	cm1_copy = deepcopy(cm1)
	cm2_copy = deepcopy(cm2)

	i = 1
	while i <= length(cm1_copy.cps)
		j = 1
		current = cm1_copy.cps[i]
		found_match = false
		while j <= length(cm2_copy.cps)
			if current.c_p.index == cm2_copy.cps[j].c_p.index &&
				current.c_e.index == cm2_copy.cps[j].c_e.index

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


	contact_point = nothing
	contact_edge = nothing

	# println(cm1_copy.cps)
	# println(cm2_copy.cps)

	if length(cm1_copy.cps) == 0

		contact_point = cm2_copy.cps[1].c_p
		contact_edge = cm2_copy.cps[1].c_e

	elseif length(cm2_copy.cps) == 0

		contact_point = cm1_copy.cps[1].c_p
		contact_edge = cm1_copy.cps[1].c_e
	else

		# println("cm1 and cm2 are not neighbors")
	end



	px = contact_point.x * cos(config.th) - contact_point.y * sin(config.th) + config.x
	py = contact_point.x * sin(config.th) + contact_point.y * cos(config.th) + config.y

	#println(px,py)

	c_p = myType.Indexed_point(contact_point.index, px, py)


	dist = dist_cp_ce(c_p, contact_edge, peg, socket)


	result1 = false
	result2 = false


	if c_p.index == id



		new_peg = Kinematics.get_peg_from_config(peg, config)
		p_c = compute_centroid(new_peg)

		px = peg.points[id].x * cos(config.th) - peg.points[id].y * sin(config.th) + config.x
		py = peg.points[id].x * sin(config.th) + peg.points[id].y * cos(config.th) + config.y


		v = [px-p_c[1], py-p_c[2]]

		F_all = [0,0]

		id_set = Array{Int}(undef, 0)

		for i = 1:length(cm1.cps)

			if id == cm1.cps[i].c_p.index
				push!(id_set, cm1.cps[i].c_e.index)
				continue
			end

			F_j = [0,0]
			F_i = deepcopy(F)

			# for j = 1:length(cm1.cps)
			# 	if j == i||id == cm1.cps[j].c_p.index
			# 		continue
			# 	end
			#
			# 	F_pj, Fnj = force_between_pair(cm1.cps[j], peg, socket, config, F)
			#
			# 	if Fnj > 0
			# 		F_j += F_pj
			# 	end
			#
			# end

			F_i = [F[1]+F_j[1], F[2]+F_j[2], F[3]]

			F_pi, Fni = force_between_pair(cm1.cps[i], peg, socket, F_i)

			if Fni > 0
				F_all += F_pi
			end

		end
		sort!(id_set)

		F_all = [F_all[1]+F[1], F_all[2]+F[2]]

		#th = compute_signed_angle([F_all[1], F_all[2], 0], [v[1], v[2], 0])

		e_id = 0

		for i = 1:length(cm2.cps)
			if id == cm2.cps[i].c_p.index
				e_id = cm2.cps[i].c_e.index

			end
		end

		vn = Array{Array{Float64}}(undef, 0)
		for i in id_set
			next = (i+1)%length(socket.points)
			if next == 0
				next = length(socket.points)
			end
			x = socket.points[i].x - socket.points[next].x
			y = socket.points[i].y - socket.points[next].y
			push!(vn, [x, y, 0])
		end

		vn[2]=-vn[2]


		th = compute_signed_angle([F_all[1], F_all[2], 0], vn[2])


		th1 = compute_signed_angle(vn[1], vn[2])/2

		p_id = Array{Int}(undef, 0)

		for i = 1:length(cm2.cps)
			push!(p_id, cm2.cps[i].c_p.index)
		end


		if (th > th1+pi/2 || th< th1-pi/2) && length(cm1.cps)==2
			result2 = true
		elseif (th > th1+pi/2 || th< th1-pi/2) && !(id in p_id)
			result1 = true
		elseif th >= th1 && e_id == id_set[1]
			result1= true
		elseif th <= th1 && e_id == id_set[2]
		    result1 = true
		end
	else
		result1 = false
		result2 = false
	end

	return result1, result2


end

function compute_centroid(peg)
	n = length(peg.points)
	A = 0
	Cx = 0
	Cy = 0
	for i = 1:n

		next = (i+1)%n
		if next == 0
			next = n
		end

		p_i = peg.points[i]
		p_n = peg.points[next]

		A += p_i.x*p_n.y-p_n.x*p_i.y


		Cx += (p_i.x+p_n.x)*(p_i.x*p_n.y-p_n.x*p_i.y )

		Cy += (p_i.y+p_n.y)*(p_i.x*p_n.y-p_n.x*p_i.y )

	end

	A = A/2
	Cx = Cx/6A
	Cy = Cy/6A

	return [Cx, Cy]

end


function is_increase_valid(F_p, cps, c_p, contact_edge)


	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, 0<=t<=1)


	p_i=contact_edge.s

	n = inward_normal(p_i, contact_edge.e)

	for i=1:length(cps)
		cp_i = cps[i].points
		ce_i = cps[i].edges


	end

	JuMP.@NLconstraint(model,
						n[1]*(c_p.x+x-p_i.x) +
						n[2]*(c_p.y+x*F_p[2]/F_p[1]-p_i.y)>=0)

	JuMP.@NLconstraint(model,
						n[1]*(c_p.x+x-p_i.x) +
						n[2]*(c_p.y+x*F_p[2]/F_p[1]-p_i.y)<=0)


	JuMP.@NLconstraint(model,
						p_i.x<=c_p.x+x<=contact_edge.e.x)


	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)

	vio = false

	p = cps.c_p
	e = cps.c_e

	if p.x+vx > e.e.x || p.x+vx < e.s.x
		vio = true
	end

	if c_p.x+vx>contact_edge.e.x || c_p.x+vx<contact_edge.s.x

		vio =true
	end

	if !all_close(n[1]*(c_p.x+vx-p_i.x) + n[2]*(c_p.y+vx*F_p[2]/F_p[1]-p_i.y), 0)

		vio =true
	end



	if status == :Optimal && !vio
		return true
	end

	return false

end


function force_between_pair(c_p::myType.contact_pair, peg::myType.myPeg,
				egs::myType.mySocketEdges, F::Array)




	e_id = c_p.c_e.index
	epi = egs.edges[e_id].s
	epj = egs.edges[e_id].e

	edge_dir = [epj.x-epi.x, epj.y-epi.y]

	i_n = inward_normal(epi, epj)

	o_n = negate_vector(i_n)

	Fn = F[1] * o_n[1] + F[2] * o_n[2]

	Fe = F[1] * edge_dir[1] + F[2] * edge_dir[2]

	if Fn <= 0
		return [F[1], F[2]], Fn
	end



	#Ft = [(CONST.friction_coeff * o_n[2] - edge_dir[2]) / (edge_dir[1] - CONST.friction_coeff * o_n[1]), 1]

	ang_e = atan(edge_dir[2], edge_dir[1])

	ang1 = ang_e - pi/2 +atan(CONST.friction_coeff)
	ang2 = ang_e - pi/2 -atan(CONST.friction_coeff)


	if F[3]<=ang1 && F[3]>=ang2

		return [0, 0], Fn

	else

		scale = abs(Fe) - abs(CONST.friction_coeff * Fn)

		F_p  = rescale_vector(edge_dir,  scale)

		if Fe > 0

			return F_p, Fn
		else
			return negate_vector(F_p), Fn
		end

	end

end

function force_between_pair(c_p::myType.contact_pair, peg::myType.myPeg,
				socket::myType.mySocket, F::Array)

	# # compute the force between any contact pairs
	# # given the orientation, force direction;
	# # and the contact point location, contact edge angle
	# # friction coefficient
	# # one can compute the combined force direction
	# # between a contact pair
	#
	# # force direction is along the given direction
	# # as the block is pushed by a constant force
	#



	e_id = c_p.c_e.index
	epi = socket.points[socket.edges[e_id].s]
	epj = socket.points[socket.edges[e_id].e]

	edge_dir = [epj.x-epi.x, epj.y-epi.y]

	i_n = inward_normal(epi, epj)

	o_n = negate_vector(i_n)

	Fn = F[1] * o_n[1] + F[2] * o_n[2]

	Fe = (F[1] * edge_dir[1] + F[2] * edge_dir[2])/vector_length(edge_dir)

	if Fn <= 0
		return [F[1], F[2]], Fn
	end



	#Ft = [(CONST.friction_coeff * o_n[2] - edge_dir[2]) / (edge_dir[1] - CONST.friction_coeff * o_n[1]), 1]

	ang_e = atan(edge_dir[2], edge_dir[1])

	ang1 = ang_e - pi/2 +atan(CONST.friction_coeff)
	ang2 = ang_e - pi/2 -atan(CONST.friction_coeff)


	if F[3]<=ang1 && F[3]>=ang2

		return [0, 0], Fn

	else

		scale = abs(Fe) - abs(CONST.friction_coeff * Fn)

		F_p  = rescale_vector(edge_dir,  scale)

		if Fe > 0

			return F_p, Fn
		else
			return negate_vector(F_p), Fn
		end

	end

end



function test()
	peg, socket = Construct_Peg_Socket.read_from_file(1)

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
	   # println(cm.cps[i].c_p.index, ", ", cm.cps[i].c_e.index)
   end

   sort_contact_mode(cm)
   # println("after sort")

   for i = 1:length(cm.cps)
	   # println(cm.cps[i].c_p.index, ", ", cm.cps[i].c_e.index)
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

    # println(cp_s[2].id, ", ", cp_s[5].id)


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
	# println(cm1_copy.cps)
	# #
	# println(cm2_copy.cps)

	# loop over cm1

end

function get_cp_cm(peg::myType.myPeg, socket::myType.mySocket,
	c_p_s::Array, model_list::Array)

	cps = Array{myType.contact_pair}(undef, 0)
	cms = Array{myType.contact_mode}(undef, 0)
	ces = myType.mySocketEdges(socket)





	for i = 1:length(c_p_s)
		cp_s = Array{myType.contact_pair}(undef, 0)
		push!(cps, myType.contact_pair(i, peg.points[c_p_s[i][1]], ces.edges[c_p_s[i][2]]))
		push!(cp_s, cps[i])
		push!(cms, myType.contact_mode(cp_s))
	end



	for i = 1:length(model_list)

		for j = 1:length(model_list[i])
			cp_s = Array{myType.contact_pair}(undef, 0)
			for k = 1:length(model_list[i][j])
				push!(cp_s, cps[model_list[i][j][k]])
			end
			push!(cms, myType.contact_mode(cp_s))
		end


	end

	return cps, cms



end


function test_vd(i)

	peg,socket = Construct_Peg_Socket.read_from_file(i)
	#peg = Adjust.get_peg_from_t(peg, socket, 0.74, 0.64)

	F = [0.,0.,0]
	f_with_p = true

	closed = false

	convex = is_convex(socket)


	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)




	config_list = cp_cf_list(peg, all_peg)

	println(c_p_s)
	println(model_list)



	Draw.draw_sth(all_peg[43])
	Draw.draw_sth(all_peg[27], (0,1,1))
	Draw.draw_open_sth(socket)


	cps, cms = get_cp_cm(peg, socket, c_p_s, model_list)


	result = false
	r1 = false
	r2 = false
	r3 = false



	if f_with_p
		F = th_get_F(config_list[43].th, 2)
	end

	if is_neighbor_modes(cms[43], cms[27])
		println("here")
		result, id = is_p_p(cms[43])
		if result
			if convex_pp(socket, cms[cmp_id[i][1]], cms[j])
				r1,r2 = p_p_direction(cms[43], cms[27], peg, socket, config_list[43], F, id)
			else
				r1,r2 = p_p_concave(cms[43], cms[27], peg, socket, config_list[43], F, id)
			end
		else
			r3 = valid_direction(cms[43], cms[27], peg, socket, config_list[43], F)
		end

	end

	println(result,r1,r2,r3)

end

function cps_trans_id(socket, F, c_p_s)


	id_set = Array{Int}(undef, 0)

	for i = 1:length(socket.edges)-1
		p_i = socket.points[socket.edges[i].s]
		n = inward_normal(p_i, socket.points[socket.edges[i].e])


		if n[1] * F[1] + n[2] * F[2] <0
			for j = 1 : length(c_p_s)
				if c_p_s[j][2] == i
					push!(id_set, j)
				end
			end

		end
	end

	return id_set


end

function get_sink(peg::myType.myPeg, socket::myType.mySocket, F::Array,
	closed = false, convex = false)

	# get sink from the insertion graph

	f_with_p = false

	if F[1]==0&&F[2]==0
		f_with_p = true
	end


	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)


	config_list = cp_cf_list(peg, all_peg)



	cps, cms = get_cp_cm(peg, socket, c_p_s, model_list)
	n = length(cms)
	g = LightGraphs.DiGraph(n+1)



	cmp_id = Array{Array}(undef, 0)

	for i = 1:n
		for j = 1:n
			result, id = is_p_p(cms[i])
			if result
				push!(cmp_id, [i,id])
				continue
			end

			if i!=j && is_neighbor_modes(cms[i], cms[j])

				if f_with_p
				 F = th_get_F(config_list[i].th, 2)
				end

				if valid_direction(cms[i], cms[j], peg, socket, config_list[i], F)

					LightGraphs.add_edge!(g, i, j)
				end
			end
		end
	end


	for i = 1:length(cmp_id)
		for j = 1: n
			if i!=j && is_neighbor_modes(cms[cmp_id[i][1]], cms[j])
				if f_with_p
					F = th_get_F(config_list[cmp_id[i][1]].th, 2)
				end


				if is_p_p_pair(cms[cmp_id[i][1]], cms[j])

					if convex_pp(socket, cms[cmp_id[i][1]], cms[j])
						result1,result2 = p_p_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F, cmp_id[i][2])
					else
						result1,result2 = p_p_concave(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F, cmp_id[i][2])
					end
					if result1

						LightGraphs.add_edge!(g, cmp_id[i][1], j)
					end

					if result2

						LightGraphs.add_edge!(g, cmp_id[i][1], n+1)

					end

				elseif valid_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F)

					LightGraphs.add_edge!(g, cmp_id[i][1], j)
				end
			end
		end

	end


	if f_with_p
		for i = 1:length(c_p_s)
			LightGraphs.add_edge!(g, n+1, i)
		end

	else
		id_set = cps_trans_id(socket, F, c_p_s)

		for i = 1:length(id_set)
			LightGraphs.add_edge!(g, n+1, i)
		end
	end



	for i =1:length(cps)
		if f_with_p
			F = th_get_F(config_list[i].th, 2)
		end
		F_p, Fn = force_between_pair(cps[i], peg, socket, F)

		if Fn<0
			LightGraphs.add_edge!(g, i, n+1)
		end
	end



	LightGraphs.collect(LightGraphs.edges(g))
	A = LightGraphs.LinAlg.adjacency_matrix(g)


	s_A = size(A, 1)
	sink = Array{Int}(undef, 0)
	for i = 1:s_A
		num = 0
		if !(i in A.rowval)
			push!(sink, i)
		end
	end


	return sink

end

function rocking_graph(peg::myType.myPeg, socket::myType.mySocket, F::Array, closed = false,
	convex = false)




	f_with_p = false

	if F[1]==0&&F[2]==0
		f_with_p = true
	end




	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)

	Draw.draw_all_peg(all_peg, socket, closed)

	config_list = cp_cf_list(peg, all_peg)

	println(c_p_s)
	println(model_list)


	cps, cms = get_cp_cm(peg, socket, c_p_s, model_list)
	n = length(cms)
	g = LightGraphs.DiGraph(n+1)



	# if valid_direction(cms[6], cms[11], peg, socket, config_list[6], F)
	# 	println("here")
	# end

	cmp_id = Array{Array}(undef, 0)

	for i = 1:n
		for j = 1:n


			if i!=j && is_neighbor_modes(cms[i], cms[j])

				LightGraphs.add_edge!(g, i, j)
			end
		end
	end




		for i = 1:length(c_p_s)
			LightGraphs.add_edge!(g, n+1, i)
		end


	# for i = 1:length(id_set)
	# 	LightGraphs.add_edge!(g, n+1, i)
	# end


	for i =1:length(cps)

			LightGraphs.add_edge!(g, i, n+1)

	end



	LightGraphs.collect(LightGraphs.edges(g))
	A = LightGraphs.LinAlg.adjacency_matrix(g)


	s_A = size(A, 1)
	sink = Array{Int}(undef, 0)
	for i = 1:s_A
		num = 0
		if !(i in A.rowval)
			push!(sink, i)
		end
	end
	println("sink:", sink)




	# if count == 0
	# 	error("no edge")
	#
	# else

	#GraphRecipes.graphplot(g, names=1:n+1,  linecolor = :darkgrey, color = :white,
	#							markersize=1.5, arrow=Plots.arrow(:closed, :head, 1, 1), shorten=0.1)

	nodesize = Vector{Float64}(undef, 0)
	nodelabelsize = Vector{Float64}(undef, 0)
	membership = Array{Int}(undef, 0)
	nodecolor = [colorant"lightseagreen", colorant"orange"]
	edgestrokec = colorant"black"

	for i = 1:n+1
		push!(nodesize, 1.)
		push!(nodelabelsize, 1.)
		push!(membership, 1)
	end

	for j in sink
		membership[j] = 2
	end
	nodefillc = nodecolor[membership]

	GraphPlot.draw(PNG("insertion_graph.png", 16cm, 16cm), GraphPlot.gplot(g, nodelabel=1:n+1,
	arrowlengthfrac=0.1,
	nodesize=nodesize, nodelabelsize=nodelabelsize, nodefillc=nodefillc,
	edgestrokec=edgestrokec, NODESIZE=0.05, NODELABELSIZE=5))
	#return g

end

function th_get_F(th,scale)
	x = scale*cos(th-pi/2)
	y = scale*sin(th-pi/2)



	return [x, y, th-pi/2]

end

function test_get_F()

	return th_get_F(0,2)

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


function model_graph(peg::myType.myPeg, socket::myType.mySocket, F::Array, closed = false,
	convex = false)

	#find insertion CMT graph


	f_with_p = false

	if F[1]==0&&F[2]==0
		f_with_p = true
	end




	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)

	Draw.draw_all_peg(all_peg, socket, closed)

	config_list = cp_cf_list(peg, all_peg)

	println(c_p_s)
	println(model_list)


	cps, cms = get_cp_cm(peg, socket, c_p_s, model_list)
	n = length(cms)
	g = LightGraphs.DiGraph(n+1)



	# if valid_direction(cms[6], cms[11], peg, socket, config_list[6], F)
	# 	println("here")
	# end

	cmp_id = Array{Array}(undef, 0)

	for i = 1:n
		for j = 1:n
			#println([i,j])
			result, id = is_p_p(cms[i])
			if result
				push!(cmp_id, [i,id])
				continue
			end

			if i!=j && is_neighbor_modes(cms[i], cms[j])
				#println("neighbors")
				# result, id = is_p_p(cms[i])
				# if result
				# 	continue
				# end
				#
				# if result
				#
				# 	result1,result2 = p_p_direction(cms[i], cms[j], peg, socket, config_list[i], F, id)
				#
				# 	if result1
				# 		count += 1
				#
				# 		LightGraphs.add_edge!(g, i, j)
				# 	end
				#
				# 	if result2
				#
				# 		LightGraphs.add_edge!(g, i, n+1)
				#
				# 	end

				if f_with_p
				 F = th_get_F(config_list[i].th, 2)
			 	end
				# println(F)
				# println(valid_direction(cms[i], cms[j], peg, socket, config_list[i], F))
				if valid_direction(cms[i], cms[j], peg, socket, config_list[i], F)

					LightGraphs.add_edge!(g, i, j)
				end
			end
		end
	end


	for i = 1:length(cmp_id)
		for j = 1: n
			if i!=j && is_neighbor_modes(cms[cmp_id[i][1]], cms[j])
				#println([cmp_id[i][1], j])
				if f_with_p
					F = th_get_F(config_list[cmp_id[i][1]].th, 2)
				end

				if is_p_p_pair(cms[cmp_id[i][1]], cms[j])

					if convex_pp(socket, cms[cmp_id[i][1]], cms[j])
						result1,result2 = p_p_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F, cmp_id[i][2])
					else
						result1,result2 = p_p_concave(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F, cmp_id[i][2])
					end

					if result1

						LightGraphs.add_edge!(g, cmp_id[i][1], j)
					end

					if result2

						LightGraphs.add_edge!(g, cmp_id[i][1], n+1)

					end
				elseif valid_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F)

					LightGraphs.add_edge!(g, cmp_id[i][1], j)
				end

			end
		end

	end


	if f_with_p
		for i = 1:length(c_p_s)
			LightGraphs.add_edge!(g, n+1, i)
		end

	else
		id_set = cps_trans_id(socket, F, c_p_s)

		for i = 1:length(id_set)
			LightGraphs.add_edge!(g, n+1, i)
		end
	end

	# for i = 1:length(id_set)
	# 	LightGraphs.add_edge!(g, n+1, i)
	# end


	for i =1:length(cps)
		if f_with_p
			F = th_get_F(config_list[i].th, 2)
		end
		F_p, Fn = force_between_pair(cps[i], peg, socket, F)

		if Fn<0
			LightGraphs.add_edge!(g, i, n+1)
		end
	end



	LightGraphs.collect(LightGraphs.edges(g))
	A = LightGraphs.LinAlg.adjacency_matrix(g)


	s_A = size(A, 1)
	sink = Array{Int}(undef, 0)
	for i = 1:s_A
		num = 0
		if !(i in A.rowval)
			push!(sink, i)
		end
	end
	println("sink:", sink)




	# if count == 0
	# 	error("no edge")
	#
	# else

	#GraphRecipes.graphplot(g, names=1:n+1,  linecolor = :darkgrey, color = :white,
	#							markersize=1.5, arrow=Plots.arrow(:closed, :head, 1, 1), shorten=0.1)

	nodesize = Vector{Float64}(undef, 0)
	nodelabelsize = Vector{Float64}(undef, 0)
	membership = Array{Int}(undef, 0)
	nodecolor = [colorant"lightseagreen", colorant"orange"]
	edgestrokec = colorant"black"

	for i = 1:n+1
		push!(nodesize, 1.)
		push!(nodelabelsize, 1.)
		push!(membership, 1)
	end

	for j in sink
		membership[j] = 2
	end
	nodefillc = nodecolor[membership]

	GraphPlot.draw(PNG("insertion_graph.png", 16cm, 16cm), GraphPlot.gplot(g, nodelabel=1:n+1,
	arrowlengthfrac=0.1,
	nodesize=nodesize, nodelabelsize=nodelabelsize, nodefillc=nodefillc,
	edgestrokec=edgestrokec, NODESIZE=0.05, NODELABELSIZE=5))
	#return g

end

function convex_pp(socket::myType.mySocket, cm1::myType.contact_mode, cm2::myType.contact_mode)
	n1 = length(cm1.cps)
	n2 = length(cm2.cps)
	egs = myType.mySocketEdges(socket)
	p_id = 0
	e_id = []

	set1 = Array{Int}(undef, 0)
	set2 = Array{Int}(undef, 0)
	for i = 1:n1
		push!(set1, cm1.cps[i].id)
	end

	for j = 1:n2
		push!(set2, cm2.cps[j].id)
	end

	set3 = setdiff(set1, set2)
	set4 = setdiff(set2, set1)

	for i = 1:n1
		if cm1.cps[i].id in set3 || cm1.cps[i].id in set4
			p_id = cm1.cps[i].c_p.index
		end
	end


	if p_id == 0
		return true
	end

	for i = 1:n1
		if cm1.cps[i].c_p.index == p_id
			push!(e_id, cm1.cps[i].c_e.index)

		end
	end

	sort!(e_id)

	e1 = egs.edges[e_id[1]]
	e2 = egs.edges[e_id[2]]

	v1 = [e1.e.x - e1.s.x, e1.e.y - e1.s.y, 0]
	v2 = [e2.e.x - e2.s.x, e2.e.y - e2.s.y, 0]

	th = compute_signed_angle(v1, v2)
	if th < 0
		return true
	else
		return false
	end

end

function th_get_F(th,scale)
	x = scale*cos(th-pi/2)
	y = scale*sin(th-pi/2)



	return [x, y, th-pi/2]

end

function test_get_F()

	return th_get_F(0,2)

end

function is_p_p(cm1)
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


function fiction_cone()

return atan(1/CONST.friction_coeff)

end

function p_p_concave(cm1::myType.contact_mode, cm2::myType.contact_mode,
						peg::myType.myPeg, socket::myType.mySocket,
						config::myType.configSE2, F::Array, id::Int)
	# point-point convex situation

	cm1_copy = deepcopy(cm1)
	cm2_copy = deepcopy(cm2)

	i = 1
	while i <= length(cm1_copy.cps)
		j = 1
		current = cm1_copy.cps[i]
		found_match = false
		while j <= length(cm2_copy.cps)
			if current.c_p.index == cm2_copy.cps[j].c_p.index &&
				current.c_e.index == cm2_copy.cps[j].c_e.index

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


	contact_point = nothing
	contact_edge = nothing

	# println(cm1_copy.cps)
	# println(cm2_copy.cps)

	if length(cm1_copy.cps) == 0

		contact_point = cm2_copy.cps[1].c_p
		contact_edge = cm2_copy.cps[1].c_e

	elseif length(cm2_copy.cps) == 0

		contact_point = cm1_copy.cps[1].c_p
		contact_edge = cm1_copy.cps[1].c_e
	else

		println("cm1 and cm2 are not neighbors")
	end



	px = contact_point.x * cos(config.th) - contact_point.y * sin(config.th) + config.x
	py = contact_point.x * sin(config.th) + contact_point.y * cos(config.th) + config.y

#println(px,py)

	c_p = myType.Indexed_point(contact_point.index, px, py)


	dist = dist_cp_ce(c_p, contact_edge, peg, socket)


	result1 = false
	result2 = false


	if c_p.index == id



		new_peg = Kinematics.get_peg_from_config(peg, config)
		p_c = compute_centroid(new_peg)

		px = peg.points[id].x * cos(config.th) - peg.points[id].y * sin(config.th) + config.x
		py = peg.points[id].x * sin(config.th) + peg.points[id].y * cos(config.th) + config.y


		v = [px-p_c[1], py-p_c[2]]

		F_all = [0,0]

		id_set = Array{Int}(undef, 0)

		for i = 1:length(cm1.cps)

			if id == cm1.cps[i].c_p.index
				push!(id_set, cm1.cps[i].c_e.index)
				continue
			end

			F_j = [0,0]
			F_i = deepcopy(F)

		# for j = 1:length(cm1.cps)
		# 	if j == i||id == cm1.cps[j].c_p.index
		# 		continue
		# 	end
		#
		# 	F_pj, Fnj = force_between_pair(cm1.cps[j], peg, socket, config, F)
		#
		# 	if Fnj > 0
		# 		F_j += F_pj
		# 	end
		#
		# end

			F_i = [F[1]+F_j[1], F[2]+F_j[2], F[3]]

			F_pi, Fni = force_between_pair(cm1.cps[i], peg, socket, F_i)

			if Fni > 0
				F_all += F_pi
			end

		end
		sort!(id_set)

		F_all = [F_all[1]+F[1], F_all[2]+F[2]]

	#th = compute_signed_angle([F_all[1], F_all[2], 0], [v[1], v[2], 0])

		e_id = 0

		for i = 1:length(cm2.cps)
			if id == cm2.cps[i].c_p.index
				e_id = cm2.cps[i].c_e.index

			end
		end

		vn = Array{Array{Float64}}(undef, 0)
		for i in id_set
			next = (i+1)%length(socket.points)
			if next == 0
				next = 4
			end
			x = socket.points[i].x - socket.points[next].x
			y = socket.points[i].y - socket.points[next].y
			push!(vn, [x, y, 0])
		end


		vn[2]=-vn[2]


		th_con = fiction_cone()

		th1 = compute_signed_angle(vn[1], [F_all[1], F_all[2], 0])

		th2= compute_signed_angle([F_all[1], F_all[2], 0], vn[2])

		p_id = Array{Int}(undef, 0)

		for i = 1:length(cm2.cps)
			push!(p_id, cm2.cps[i].c_p.index)
		end


		if th1 < 0 && th2 < 0 && length(cm1.cps)==2
			result2 = true
		elseif th1 < 0 && th2 < 0  && !(id in p_id)
			result1 = true
		elseif th1 < th_con && e_id == id_set[1]
			result1= true
		elseif th2 < th_con && e_id == id_set[2]
			result1 = true
		end



	else
		result1 = false
		result2 = false
	end

	return result1, result2
end


function is_p_p_pair(cm1, cm2)
	# test if the different CP can cause point-point situation

	n = length(cm1.cps)
	m = length(cm2.cps)

	set1 = Array{Int}(undef, 0)
	set2 = Array{Int}(undef, 0)
	for i = 1:n
		push!(set1, cm1.cps[i].c_p.index)
	end

	for j = 1:m
		push!(set2, cm2.cps[j].c_p.index)
	end

	set3 = setdiff(set1, set2)
	set4 = setdiff(set2, set1)

	if length(set3) == 0 && length(set4) == 0
		return true
	end

	return false


end


function compute_config(peg1::myType.myPeg, peg2::myType.myPeg)
	#find the configuration which transfer peg 1 to peg 2


	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)

	for i = 1:2


		p1 = peg1.points[i]
		p2 = peg2.points[i]



		JuMP.@NLconstraint(model,
			(cos(th) * p1.x - sin(th) * p1.y + x)  >= p2.x);

		JuMP.@NLconstraint(model,
			(cos(th) * p1.x - sin(th) * p1.y + x)  <= p2.x);

		JuMP.@NLconstraint(model,
			(sin(th) * p1.x + cos(th) * p1.y + y)  >= p2.y);

		JuMP.@NLconstraint(model,
			(sin(th) * p1.x + cos(th) * p1.y + y)  <= p2.y);


	end

	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	config = configSE2(vx,vy,vth)
	new_peg = Kinematics.get_peg_from_config(peg1,config)

	# Draw.draw_sth(peg1)
	# Draw.draw_sth(new_peg,(0,1,1),2.0)
	# Draw.draw_sth(peg2,(1,1,0),2.0)

	vio = false






		# for j = 1:length(peg1.points)
		#
		# 	p1 = peg1.points[j]
		# 	p2 = peg2.points[j]
		#
		#
		#
		# 	# println(n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
		# 	# 	n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) )
		# 	if !all_close((cos(vth) * p1.x - sin(vth) * p1.y + vx), p2.x) ||
		# 		!all_close((sin(vth) * p1.x + cos(vth) * p1.y + vy), p2.y)
		#
		# 		vio = true
		# 	end
		#
		# end






	if status == :Optimal && !vio
		return vx, vy, vth
	end

	return nothing, nothing, nothing





end


function cp_cf_list(peg1::myType.myPeg, pl::Array{myType.myPeg})
	# find confiuration list for all CP


	config_list = Array{myType.configSE2}(undef, 0)

	for i = 1:length(pl)
		vx,vy,vth = compute_config(peg1, pl[i])

		if vx == nothing
			error("no config")
		end
		push!(config_list, myType.configSE2(vx,vy,vth))
	end

	return config_list


end

function is_convex(socket::myType.mySocket)
	#test if the socket is convex

	n1 = length(socket.points)

	ps = Array{Array{Float64,1}}(undef, 0)
	for i = 1:n1
		push!(ps, [socket.points[i].x, socket.points[i].y])
	end

	hull_ps = Convexhull.convex_hull(ps)



	n2 = length(hull_ps)


	if n1 == n2
		return true
	end

	return false

end


function get_graph_edges(peg::myType.myPeg, socket::myType.mySocket, F=[0,0,0],
	closed = false, convex = false)

	# return the CMT graph


	f_with_p = false

	if F[1]==0&&F[2]==0
		f_with_p = true
	end



	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)


	config_list = cp_cf_list(peg, all_peg)



	cps, cms = get_cp_cm(peg, socket, c_p_s, model_list)
	n = length(cms)
	g = LightGraphs.DiGraph(n+1)


	cmp_id = Array{Array}(undef, 0)

	for i = 1:n
		for j = 1:n

			result, id = is_p_p(cms[i])
			if result
				push!(cmp_id, [i,id])
				continue
			end

			if i!=j && is_neighbor_modes(cms[i], cms[j])


				if f_with_p
				 F = th_get_F(config_list[i].th, 2)
			 	end

				if valid_direction(cms[i], cms[j], peg, socket, config_list[i], F)

					LightGraphs.add_edge!(g, i, j)
				end
			end
		end
	end


	for i = 1:length(cmp_id)
		for j = 1: n
			if i!=j && is_neighbor_modes(cms[cmp_id[i][1]], cms[j])

				if f_with_p
					F = th_get_F(config_list[cmp_id[i][1]].th, 2)
				end

				if is_p_p_pair(cms[cmp_id[i][1]], cms[j])

					if convex_pp(socket, cms[cmp_id[i][1]], cms[j])
						result1,result2 = p_p_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F, cmp_id[i][2])
					else
						result1,result2 = p_p_concave(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F, cmp_id[i][2])
					end

					if result1

						LightGraphs.add_edge!(g, cmp_id[i][1], j)
					end

					if result2

						LightGraphs.add_edge!(g, cmp_id[i][1], n+1)

					end
				elseif valid_direction(cms[cmp_id[i][1]], cms[j], peg, socket, config_list[cmp_id[i][1]], F)

					LightGraphs.add_edge!(g, cmp_id[i][1], j)
				end

			end
		end

	end


	if f_with_p
		for i = 1:length(c_p_s)
			LightGraphs.add_edge!(g, n+1, i)
		end

	else
		id_set = cps_trans_id(socket, F, c_p_s)

		for i = 1:length(id_set)
			LightGraphs.add_edge!(g, n+1, i)
		end
	end




	for i =1:length(cps)
		if f_with_p
			F = th_get_F(config_list[i].th, 2)
		end
		F_p, Fn = force_between_pair(cps[i], peg, socket, F)

		if Fn<0
			LightGraphs.add_edge!(g, i, n+1)
		end
	end



	g = LightGraphs.collect(LightGraphs.edges(g))


	return g

end

function test_graph(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)

	F = [0,0,0]


	closed = false

	convex = is_convex(socket)


	sink = get_sink(peg,socket, F, closed)
	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)


	#peg = Adjust.get_peg_from_t(peg, socket, 0.74, 0.64)
	# for i = 1:length(sink)
	# 	PyPlot.figure()
	# 	Draw.draw_open_sth(socket)
	# 	Draw.draw_sth(peg,[1,0,1],2.0)
	# 	Draw.draw_sth(all_peg[sink[i]],[0,1,1],2.0)
	# end

	model_graph(peg,socket, F, closed, convex)

end

function test_rocking_graph(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)

	p = compute_centroid(peg)


	F = [0,0,0]


	closed = true

	convex = is_convex(socket)



	all_peg, c_p_s, model_list = all_combinations(peg, socket, closed, convex)


	#peg = Adjust.get_peg_from_t(peg, socket, 0.74, 0.64)
	# for i = 1:length(sink)
	# 	PyPlot.figure()
	# 	Draw.draw_open_sth(socket)
	# 	Draw.draw_sth(peg,[1,0,1],2.0)
	# 	Draw.draw_sth(all_peg[sink[i]],[0,1,1],2.0)
	# end

	rocking_graph(peg,socket, F, closed, convex)

end





end
