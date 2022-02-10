module Rotation
import LinearAlgebra
import JuMP

import MathOptInterface
import NLopt


import ..CONST
import ..Kinematics
import ..Constraints
import ..Construct_Peg_Socket
import ..Insertion
import ..Convexhull
import ..Partial_order
import ..Draw

import PyPlot

#import ..myType
using ..myType
using LazySets
using Plots

############### this module is to rotate peg by move t ########


function all_close(a, b)
	if abs(a - b) < 0.011
		return true
	end
	return false;

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
        # println("here")
		return nothing
	end


	# find solution
	t1 = (-B + sqrt(B^2-4*A*C)) / (2 * A)
	t2 = (-B - sqrt(B^2-4*A*C)) / (2 * A)

    if t1 <= 1 && t1>=0
        return t1
    elseif t2<=1 && t2>=0
        return t2
    end

	return nothing

end


function get_t(e::Indexed_edge, p::Indexed_point)
    return (p.x-e.s.x)/(e.e.x-e.s.x)
end

function get_xy(e::Indexed_edge, t::Float64)
    v = [e.e.x-e.s.x, e.e.y-e.s.y]
    p = [e.s.x, e.s.y]
    p = p+v*t
    return myPoint2d(p[1], p[2])
end

function get_xy(e::Indexed_edge, p::Indexed_point, t::Float64)
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

function pockets_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(-pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function signed_distance(e::myType.Indexed_edge, p::myType.Indexed_point)
    v = [p.x-e.s.x, p.y-e.e.y]
    n = inward_normal(e.s, e.e)
    return v[1]*n[1]+v[2]*n[2]
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

function vector_length(v::Array)
	scale = 0
	for i = 1:length(v)
		scale += v[i]*v[i]
	end

	return sqrt(scale)

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

function uncontact_dis(es::myType.mySocketEdges, p_set::Array{Indexed_point}, cp_id::Array{Int})
    vio_dis = Array{Array{Float64}}(undef, 0)
    ne = length(es.edges)
    np = length(p_set)
    for i = 1:np
        if p_set[i].index in cp_id
            continue
        end
        dis = Array{Float64}(undef, 0)
        for j = 1:ne
            d = signed_distance(es.edges[j], p_set[i])
            push!(dis, d)
        end
        push!(vio_dis, dis)
    end
    return vio_dis

end


function test_vio(peg::myType.myPeg, socket::myType.mySocket)

	convex = Insertion.is_convex(socket)
	vio = false

	if !convex
	    hull, pockets = Convexhull.get_convexpocket(socket)



		vio_id = Array{Array{Int}}(undef, 0)

	    for i = 1:length(hull.edges)
	        p_i = socket.points[hull.edges[i].s]
	        # get the inward normal
	        n = inward_normal(p_i, socket.points[hull.edges[i].e])

	        for j = 1:length(peg.points)

	            # for each potential contact point on peg
	            # add constraint of the positive
	            p_j = peg.points[j]
	            vn = [p_j.x - p_i.x, p_j.y - p_i.y]
	             # println(n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
	             # 	n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) )
	            if vn[1]*n[1]+vn[2]*n[2] < -0.0001

	                vio = true
	            end

	        end
	    end

		for i=1:length(pockets)

			for j = 1:length(peg.points)


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

				if dis_set[1][1]<= 0.0001
					continue
				end

				vio = true



				if !([dis_set[1][2],dis_set[1][3]] in vio_id)
					push!(vio_id, [dis_set[1][2],dis_set[1][3]])
				end

			end
		end
	else
		for i = 1:length(socket.edges)
	        p_i = socket.points[socket.edges[i].s]
	        # get the inward normal
	        n = inward_normal(p_i, socket.points[socket.edges[i].e])

	        for j = 1:length(peg.points)

	            # for each potential contact point on peg
	            # add constraint of the positive
	            p_j = peg.points[j]
	            vn = [p_j.x - p_i.x, p_j.y - p_i.y]
	             # println(n[1] * ((cos(vth) * cp_j.x - sin(vth) * cp_j.y + vx) - p_i.x) +
	             # 	n[2] * ((sin(vth) * cp_j.x + cos(vth) * cp_j.y + vy) - p_i.y) )
	            if vn[1]*n[1]+vn[2]*n[2] < -0.0001

	                vio = true
	            end

	        end
	    end
	end

	return vio


end

function compute_config(p1::Array{myPoint2d}, p2::Array{Indexed_point})

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)

    n = length(p1)

    for i = 1:n

	       JuMP.@NLconstraint(model,
		         (cos(th) * p2[i].x - sin(th) * p2[i].y + x)  >= p1[i].x);

	       JuMP.@NLconstraint(model,
		         (cos(th) * p2[i].x - sin(th) * p2[i].y + x)  <= p1[i].x);

	       JuMP.@NLconstraint(model,
		         (sin(th) * p2[i].x + cos(th) * p2[i].y + y)  >= p1[i].y);

           JuMP.@NLconstraint(model,
      		     (sin(th) * p2[i].x + cos(th) * p2[i].y + y)  <= p1[i].y);

    end





	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	config = configSE2(vx,vy,vth)

	# new_peg = Kinematics.get_peg_from_config(peg1,config)
    #
    #
	# Draw.draw_sth(new_peg,(0,1,1),2.0)


	vio = false







    for i = 1:n

        # println([(cos(vth) * p2[i].x - sin(vth) * p2[i].y + vx), (sin(vth) * p2[i].x + cos(vth) * p2[i].y + vy)])
        # println([p1[i].x, p1[i].y])
	       if !all_close((cos(vth) * p2[i].x - sin(vth) * p2[i].y + vx), p1[i].x) ||
		      !all_close((sin(vth) * p2[i].x + cos(vth) * p2[i].y + vy), p1[i].y)
			  # println([i, abs((cos(vth) * p2[i].x - sin(vth) * p2[i].y + vx)-p1[i].x)
			  # ,abs((sin(vth) * p2[i].x + cos(vth) * p2[i].y + vy)-p1[i].y)])
		      vio = true
	       end
    end







	if !vio
		return config
	end


    # println("no config")
	return nothing





end


function compute_config(peg::myType.myPeg, p::Array{myPoint2d}, id::Array)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)

    n = length(id)

    for i = 1:n

        p_i = peg.points[id[i]]

	       JuMP.@NLconstraint(model,
		         (cos(th) * p_i.x - sin(th) * p_i.y + x)  >= p[i].x);

	       JuMP.@NLconstraint(model,
		         (cos(th) * p_i.x - sin(th) * p_i.y + x)  <= p[i].x);

	       JuMP.@NLconstraint(model,
		         (sin(th) * p_i.x + cos(th) * p_i.y + y)  >= p[i].y);

           JuMP.@NLconstraint(model,
      		     (sin(th) * p_i.x + cos(th) * p_i.y + y)  <= p[i].y);

    end





	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	config = configSE2(vx,vy,vth)

	# new_peg = Kinematics.get_peg_from_config(peg1,config)
    #
    #
	# Draw.draw_sth(new_peg,(0,1,1),2.0)


	vio = false







    for i = 1:n
        p_i = peg.points[id[i]]
        # println([(cos(vth) * p2[i].x - sin(vth) * p2[i].y + vx), (sin(vth) * p2[i].x + cos(vth) * p2[i].y + vy)])
        # println([p1[i].x, p1[i].y])
	       if !all_close((cos(vth) * p_i.x - sin(vth) *p_i.y + vx), p[i].x) ||
		      !all_close((sin(vth) * p_i.x + cos(vth) * p_i.y + vy), p[i].y)
			  # println([i, abs((cos(vth) * p2[i].x - sin(vth) * p2[i].y + vx)-p1[i].x)
			  # ,abs((sin(vth) * p2[i].x + cos(vth) * p2[i].y + vy)-p1[i].y)])
		      vio = true
	       end
    end







	if !vio
		return config
	end


    #println("no config")
	return nothing





end

function compute_config(p1::Array, p2::Array, ip1::Indexed_point, ip2::Indexed_point)

	model = JuMP.Model(solver=NLopt.NLoptSolver(algorithm=:LN_COBYLA, constrtol_abs=0.00001))

	JuMP.@variable(model, x)
	JuMP.@variable(model, y)
	JuMP.@variable(model, th)


	JuMP.@NLconstraint(model,
		(cos(th) * ip1.x - sin(th) * ip1.y + x)  >= p1[1]);

	JuMP.@NLconstraint(model,
		(cos(th) * ip1.x - sin(th) * ip1.y + x)  <= p1[1]);

	JuMP.@NLconstraint(model,
		(sin(th) * ip1.x + cos(th) * ip1.y + y)  >= p1[2]);

	JuMP.@NLconstraint(model,
		(sin(th) * ip1.x + cos(th) * ip1.y + y)  <= p1[2]);


    JuMP.@NLconstraint(model,
    	(cos(th) * ip2.x - sin(th) * ip2.y + x)  >= p2[1]);

	JuMP.@NLconstraint(model,
    	(cos(th) * ip2.x - sin(th) * ip2.y + x)  <= p2[1]);

	JuMP.@NLconstraint(model,
    	(sin(th) * ip2.x + cos(th) * ip2.y + y)  >= p2[2]);

	JuMP.@NLconstraint(model,
		(sin(th) * ip2.x + cos(th) * ip2.y + y)  <= p2[2]);




	JuMP.@NLobjective(model, Min, 1==1)


	status = JuMP.solve(model)


	vx = JuMP.getvalue(x)
	vy = JuMP.getvalue(y)
	vth = JuMP.getvalue(th)

	config = configSE2(vx,vy,vth)
	# new_peg = Kinematics.get_peg_from_config(peg1,config)

	# Draw.draw_sth(peg1)
	# Draw.draw_sth(new_peg,(0,1,1),2.0)
	# Draw.draw_sth(peg2,(1,1,0),2.0)

	vio = false








	if !all_close((cos(vth) * ip1.x - sin(vth) * ip1.y + vx), p1[1]) ||
		!all_close((sin(vth) * ip1.x + cos(vth) * ip1.y + vy), p1[2])

		vio = true
	end

    if !all_close((cos(vth) * ip2.x - sin(vth) * ip2.y + vx), p2[1]) ||
        !all_close((sin(vth) * ip2.x + cos(vth) * ip2.y + vy), p2[2])

        vio = true
    end







	if status == :Optimal && !vio
		return config
	end

    println("no config")
	return nothing





end


function rotate_edge(socket::myType.mySocket, th::Float64, e_id::Int, mode::Int)
	n = length(socket.points) - 1

	# line_list = Array{Array{Array}}(undef, 0)
	#
	# for i = 1:n-1
	#     push!(line_list, [[socket.points[i].x,socket.points[i].y], [socket.points[i+1].x, socket.points[i+1].y]])
	# end



	if e_id <= n/2


		p_s = socket.points[socket.edges[e_id].s]
		p_e = socket.points[socket.edges[e_id].e]




		R = Kinematics.rotation_matrix_2D(th)

		ps = [p_s.x, p_s.y]
		pe = [p_e.x, p_e.y]

		next = (e_id+1)%n
		last = (e_id-1)%n
		if next == 0
			next = n
		end
		if last == 0
			last = n
		end

		pn = [socket.points[socket.edges[next].e].x, socket.points[socket.edges[next].e].y]

		if mode == 0
			pc = ps
			pns = [socket.points[socket.edges[next].s].x, socket.points[socket.edges[next].s].y]
			pne = [socket.points[socket.edges[next].e].x, socket.points[socket.edges[next].e].y]
		else mode == 1
			pc = pe
			pns = [socket.points[socket.edges[last].s].x, socket.points[socket.edges[last].s].y]
			pne = [socket.points[socket.edges[last].e].x, socket.points[socket.edges[last].e].y]
		end

		pe_new = pc + R*(pe-pc)
		ps_new = pc + R*(ps-pc)



		l1 = get_line(ps_new,pe_new)
		l2 = get_line(pns,pne)

		new_socket = deepcopy(socket)

		if mode == 0
			pe_new = get_intersection(l1,l2)
			new_socket.points[socket.edges[e_id].e].x = pe_new[1]
			new_socket.points[socket.edges[e_id].e].y = pe_new[2]
		else
			ps_new = get_intersection(l1,l2)
			new_socket.points[socket.edges[e_id].s].x = ps_new[1]
			new_socket.points[socket.edges[e_id].s].y = ps_new[2]
		end




	# vn = pe_new - pe


	# for i = e_id+1:n
	# 	new_socket.points[i].x = socket.points[i].x + vn[1]
	# 	new_socket.points[i].y = socket.points[i].y + vn[2]
	# end
	else

		p_s = socket.points[socket.edges[e_id].s]
		p_e = socket.points[socket.edges[e_id].e]




		R = Kinematics.rotation_matrix_2D(th)

		ps = [p_s.x, p_s.y]
		pe = [p_e.x, p_e.y]

		next = (e_id+1)%n
		last = (e_id-1)%n
		if next == 0
			next = n

		end
		if last == 0
			last = n
		end



		if mode == 0
			pc = pe
			pns = [socket.points[socket.edges[last].s].x, socket.points[socket.edges[last].s].y]
			pne = [socket.points[socket.edges[last].e].x, socket.points[socket.edges[last].e].y]
		else mode == 1
			pc = ps
			pns = [socket.points[socket.edges[next].s].x, socket.points[socket.edges[next].s].y]
			pne = [socket.points[socket.edges[next].e].x, socket.points[socket.edges[next].e].y]
		end


		pe_new = pc + R*(pe-pc)
		ps_new = pc + R*(ps-pc)



		l1 = get_line(ps_new,pe_new)
		l2 = get_line(pns,pne)



		new_socket = deepcopy(socket)

		if mode == 0
			ps_new = get_intersection(l1,l2)
			new_socket.points[socket.edges[e_id].s].x = ps_new[1]
			new_socket.points[socket.edges[e_id].s].y = ps_new[2]
		else
			pe_new = get_intersection(l1,l2)
			new_socket.points[socket.edges[e_id].e].x = pe_new[1]
			new_socket.points[socket.edges[e_id].e].y = pe_new[2]
		end

	end

	return new_socket

end


function rotate_socket(socket::myType.mySocket, r::Array, i::Int, th::Float64)
	n_s = length(socket.points) - 1

	mode = r[i]

	socket = rotate_edge(socket, th, i, mode)
	socket = rotate_edge(socket, -th, n_s + 1 -i, mode)

	return socket

end

function capture_rotate(socket::myType.mySocket, th::Float64)
	n_s = length(socket.points) - 1


	socket = rotate_edge(socket, th, 1, 0)
	socket = rotate_edge(socket, -th, n_s  , 0)

	return socket

end

function test_rotatesocket(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	new_socket = adj_upsocket(socket)
	#Draw.draw_sth(peg)
	Draw.draw_open_sth(socket)
	Draw.draw_open_sth(new_socket)


end

function adj_upsocket(socket::myType.mySocket, th=0.1)
    n = length(socket.points) - 1

	socket = rotate_edge(socket, th, 1, 0)
	socket = rotate_edge(socket, -th, n, 0)
	return socket
end

function adj_lowsocket(socket::myType.mySocket, th=0.01)
    n = length(socket.points) - 1

    if n <= 4
        error("no lowedge")
    end


	socket = rotate_edge(socket, -th, 2)
	socket = rotate_edge(socket, th, n-2)
	return socket

end

function MoC_CPindex(cm::myType.contact_mode)
	n = length(cm.cps)
	cp_id = []
	for i = 1:n
		pid = cm.cps[i].c_p.index
		eid = cm.cps[i].c_e.index
		push!(cp_id, [pid, eid])

	end

	return cp_id


end

function get_rotate_th(peg::myType.myPeg, socket::myType.mySocket,
	F = [0,0,0], closed = false, convex = false, pic=false)




    all_peg, c_p_s, model_list = Insertion.all_combinations(peg, socket, closed, convex)

    cps, cms = Insertion.get_cp_cm(peg, socket, c_p_s, model_list)

    config_list =Insertion.cp_cf_list(peg, all_peg)

    maxth_set = []

    poset = Partial_order.get_poset(peg, socket, closed, convex)

	if poset == nothing
		return nothing, nothing, nothing
	end
	# println(poset)
    for i in poset
        max_th = rotate_th(socket, peg, cms[i], config_list[i], pic)
        #println([i, th1, th2])
		# println(max_th)
		cp_id = MoC_CPindex(cms[i])

        push!(maxth_set, [max_th, cp_id])
    end

    ths = deepcopy(maxth_set)

    sort!(maxth_set, by = x -> x[1])

    return ths, maxth_set[length(maxth_set)], maxth_set[length(maxth_set)-1]
end


function rotate_th(socket::myType.mySocket, peg::myType.myPeg,
    cm::myType.contact_mode, config::configSE2, pic=false)
    n = length(cm.cps)
    th = Array{Float64}(undef, 0)

    id_set = Array{Array{Int}}(undef, 0)

    th1_set = Array{Float64}(undef, 0)
    th2_set = Array{Float64}(undef, 0)
    max_th = 0


	if pic
     	PyPlot.figure()
		PyPlot.axis("off")
	end


    th1 = rotate_left(socket, peg, cm, config, pic)

    th2 = rotate_right(socket, peg, cm, config, pic)
    new_peg = Kinematics.get_peg_from_config(peg, config)

	if pic
        Draw.draw_sth(socket, (0,0,0), 3.)
        # Draw.draw_sth(new_peg, (1,1,0), 2.)
        # Draw.draw_sth(peg, (1,0,1), 2.)
	end

    max_th = max(abs(th1), abs(th2))
    return max_th

end


function rotate_left(socket::myType.mySocket, peg::myType.myPeg, cm::myType.contact_mode,
	config::myType.configSE2, pic=false)
    new_peg = Kinematics.get_peg_from_config(peg, config)


    np = length(new_peg.points)
    ne = length(socket.edges)

    p_set = Array{Indexed_point}(undef, 0)
    for i = 1:np
        push!(p_set, new_peg.points[i])
    end

    es = myType.mySocketEdges(socket)

    cp_id = Array{Int}(undef, 0)
    ce_id = Array{Int}(undef, 0)
    cp_set = Array{Indexed_point}(undef, 0)
    ce_set = Array{Indexed_edge}(undef, 0)
    for i = 1:length(cm.cps)
        push!(cp_id, cm.cps[i].c_p.index)
        push!(ce_id, cm.cps[i].c_e.index)
        push!(cp_set, new_peg.points[cm.cps[i].c_p.index])
        push!(ce_set, es.edges[cm.cps[i].c_e.index])
    end

    ncp = length(cp_set)
    nce = length(ce_set)

    vio_dis = uncontact_dis(es, p_set, cp_id)

    if np == 1
        return nothing
    end

    rs = Array{Float64}(undef, 0)
    for i = 1:ncp-1
        push!(rs, distance(cp_set[1], cp_set[i+1]))
    end


    t1 = get_t(ce_set[1], cp_set[1])
    vio = false
    th = 0
    new_config=myType.configSE2(0,0,0)
    config1=myType.configSE2(0,0,0)
    while t1<0.99
        new_peg = Kinematics.get_peg_from_config(peg, config)

        p = get_xy(ce_set[1], t1)

        newcp_set = Array{myPoint2d}(undef, 0)

        push!(newcp_set, p)
        t = Array{Float64}(undef, 0)
        for j = 1:length(rs)
            tj = line_circle_intersection(ce_set[j+1].s, ce_set[j+1].e, p, rs[j])
            if tj == nothing
                if pic
                    Draw.draw_sth(new_peg, (0,0,1), 3.)
                end
                return config.th

            end
            push!(t, tj)
            pj = get_xy(ce_set[j+1], t[j])
            push!(newcp_set, pj)
        end



        new_config = compute_config(peg, newcp_set, cp_id)

        if new_config == nothing
            if pic
                Draw.draw_sth(new_peg, (0,0,1), 3.)
            end

            return config.th
        else
            config1 = new_config
        end



        peg2 = Kinematics.get_peg_from_config(peg, config1)

        pj_set = Kinematics.compute_points_from_config(peg, config1)

		# if pic
		# 	Draw.draw_sth(peg2)
		# end


        for k = 1:np
            if test_vio(peg2, socket)
                if pic
                    Draw.draw_sth(new_peg, (0,0,1), 3.)
                end
                return config.th
            end
        end



        t1 += 0.01
        #println(t1)
        config = config1
    end

    peg3 = Kinematics.get_peg_from_config(peg, config)
	if pic
		Draw.draw_sth(peg3, (0,0,1), 3.)
	end
    return config.th




end

function rotate_right(socket::myType.mySocket, peg::myType.myPeg, cm::myType.contact_mode,
	config::myType.configSE2, pic=false)
    new_peg = Kinematics.get_peg_from_config(peg, config)



        np = length(new_peg.points)
        ne = length(socket.edges)

        p_set = Array{Indexed_point}(undef, 0)
        for i = 1:np
            push!(p_set, new_peg.points[i])
        end

        es = myType.mySocketEdges(socket)

        cp_id = Array{Int}(undef, 0)
        ce_id = Array{Int}(undef, 0)
        cp_set = Array{Indexed_point}(undef, 0)
        ce_set = Array{Indexed_edge}(undef, 0)
        for i = 1:length(cm.cps)
            push!(cp_id, cm.cps[i].c_p.index)
            push!(ce_id, cm.cps[i].c_e.index)
            push!(cp_set, new_peg.points[cm.cps[i].c_p.index])
            push!(ce_set, es.edges[cm.cps[i].c_e.index])
        end

        ncp = length(cp_set)
        nce = length(ce_set)

        vio_dis = uncontact_dis(es, p_set, cp_id)

        if np == 1
            return nothing
        end

        rs = Array{Float64}(undef, 0)
        for i = 1:ncp-1
            push!(rs, distance(cp_set[1], cp_set[i+1]))
        end


        t1 = get_t(ce_set[1], cp_set[1])
        vio = false
        th = 0
        new_config=myType.configSE2(0,0,0)
        config1=myType.configSE2(0,0,0)
        while t1>0
            new_peg = Kinematics.get_peg_from_config(peg, config)

            p = get_xy(ce_set[1], t1)

            newcp_set = Array{myPoint2d}(undef, 0)

            push!(newcp_set, p)
            t = Array{Float64}(undef, 0)
            for j = 1:length(rs)
                tj = line_circle_intersection(ce_set[j+1].s, ce_set[j+1].e, p, rs[j])
                if tj == nothing
                    if pic
                        Draw.draw_sth(new_peg, (0,0,1), 3.)
                    end
                    return config.th

                end
                push!(t, tj)
                pj = get_xy(ce_set[j+1], t[j])
                push!(newcp_set, pj)
            end



            new_config = compute_config(peg, newcp_set, cp_id)

            if new_config == nothing
                if pic
                    Draw.draw_sth(new_peg, (0,0,1), 3.)
                end

                return config.th
            else
                config1 = new_config
            end



            peg2 = Kinematics.get_peg_from_config(peg, config1)

            pj_set = Kinematics.compute_points_from_config(peg, config1)

    		# if pic
    		# 	Draw.draw_sth(peg2)
    		# end


            for k = 1:np
                if test_vio(peg2, socket)
                    if pic
                        Draw.draw_sth(new_peg, (0,0,1), 3.)
                    end
                    return config.th
                end
            end



            t1 += 0.01
            #println(t1)
            config = config1
        end

        peg3 = Kinematics.get_peg_from_config(peg, config)
    	if pic
    		Draw.draw_sth(peg3, (0,0,1), 3.)
    	end
        return config.th



end




function test_rotate(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
    F = [0,0,0]
    closed = true
    convex = Insertion.is_convex(socket)
    pic = true
	thset, fth, sth = get_rotate_th(peg, socket, F, closed, convex, pic)
	println(thset)
	println(fth)
	println(sth)
end














end
