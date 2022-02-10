module Opjoint


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
import ..Adjust
import ..Rotation
import ..Stability_gradient
import ..Insertion_gradient


import PyPlot
import Plots


function compute_signed_angle(v1, v2)

	vn = [0,0,1]

	#os = LinearAlgebra.dot(v1,v2)/v1[1]^2+v2

	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return th


end

function compute_abs_angle(v1, v2)

	vn = [0,0,1]

	#os = LinearAlgebra.dot(v1,v2)/v1[1]^2+v2

	tan = LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn)/LinearAlgebra.dot(v1,v2)


	th = atan(LinearAlgebra.dot(LinearAlgebra.cross(v1,v2), vn), LinearAlgebra.dot(v1,v2))



	return abs(th)


end

function pockets_normal(point1::Indexed_point, point2::Indexed_point)

	# assume the points are given counter-clockwise

	R = Kinematics.rotation_matrix_2D(-pi/2)

	v = [point2.x-point1.x, point2.y-point1.y]

	result = R * v

	return LinearAlgebra.normalize(result)

end

function all_close(a, b)
	if abs(a - b) < 0.011
		return true
	end
	return false;

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
    return [A, B, C]
end

function get_line(p1::Indexed_point, p2::Indexed_point)
    A = (p1.y - p2.y)
    B = (p2.x - p1.x)
    C = (p1.x*p2.y - p2.x*p1.y)
    return [A, B, C]
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

function get_xy(e::Indexed_edge, p::Indexed_point, ::Float64)
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

function vio_test(es::myType.mySocketEdges, p1::Indexed_point, p2::Indexed_point)

    ne = length(es.edges)

    for i = 1:ne

        d1 = signed_distance(es.edges[i], p1)
        d2 = signed_distance(es.edges[i], p2)

        if d1*d2 <= 0
            return true
        end
    end

    return false

end

function test_vio(peg::myType.myPeg, socket::myType.mySocket)

    hull, pockets = Convexhull.get_convexpocket(socket)

	vio = false

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

	return vio


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


    println("no config")
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
    while t1<0.99 && vio == false

        p = get_xy(ce_set[1], t1)

        newcp_set = Array{myPoint2d}(undef, 0)

        push!(newcp_set, p)
        t = Array{Float64}(undef, 0)
        for j = 1:length(rs)
            tj = line_circle_intersection(ce_set[j+1].s, ce_set[j+1].e, p, rs[j])
            if tj == nothing

                vio = true
                break
            end
            push!(t, tj)
            pj = get_xy(ce_set[j+1], t[j])
            push!(newcp_set, pj)
        end

        if vio == true
            break
        end

        config = compute_config([newcp_set[1], newcp_set[2]], [cp_set[1], cp_set[2]])
        peg2 = Kinematics.get_peg_from_config(new_peg, config)

        pj_set = Kinematics.compute_points_from_config(new_peg, config)

		# if pic
		# 	Draw.draw_sth(peg2)
		# end


        for k = 1:np
            if test_vio(peg2, socket)
                vio = true
                break
            end
        end

        if vio == true
            break
        end


        t1 += 0.01
        #println(t1)
        if vio == false
            new_config = config


        end

    end

    peg3 = Kinematics.get_peg_from_config(new_peg, new_config)
	if pic
		Draw.draw_sth(peg3, (0,1,0), 3.)
	end
    return new_config.th




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
    while t1>0 && vio == false

        p = get_xy(ce_set[1], t1)

        newcp_set = Array{myPoint2d}(undef, 0)
        push!(newcp_set, p)
        t = Array{Float64}(undef, 0)
        for j = 1:length(rs)
            tj = line_circle_intersection(ce_set[j+1].s, ce_set[j+1].e, p, rs[j])
            if tj == nothing

                vio = true

                break
            end
            push!(t, tj)
            pj = get_xy(ce_set[j+1], t[j])
            push!(newcp_set, pj)
        end

        if vio == true
            break
        end

        config = compute_config([newcp_set[1], newcp_set[2]], [cp_set[1], cp_set[2]])
        peg2 = Kinematics.get_peg_from_config(new_peg, config)

		# if pic
		# 	Draw.draw_sth(peg2)
		# end
        pj_set = Kinematics.compute_points_from_config(new_peg, config)

        for k = 1:np

            if test_vio(peg2, socket)
                vio = true

                break
            end
        end

        if vio == true
            break
        end


        t1 -= 0.01
        #println(t1)
        if vio == false
            new_config = config


        end

    end

    peg3 = Kinematics.get_peg_from_config(new_peg, new_config)
	if pic
     	Draw.draw_sth(peg3,(0,0,1),3.)
 	end
    return new_config.th




end

function rotate_th(socket::myType.mySocket, peg::myType.myPeg, cm::myType.contact_mode,
	config::configSE2, pic=false)
    n = length(cm.cps)
    th = Array{Float64}(undef, 0)

    id_set = Array{Array{Int}}(undef, 0)

    th1_set = Array{Float64}(undef, 0)
    th2_set = Array{Float64}(undef, 0)
    max_th = 0
    for i = 1:n

        for j = 1:n
            if j == i
                continue
            end

            v = sort!([i, j])
            if v in id_set
                continue
            end

			if pic
             	PyPlot.figure()
				PyPlot.axis("off")
			end
            newcm = myType.contact_mode([cm.cps[i], cm.cps[j]])

            th1 = rotate_left(socket, peg, newcm, config, pic)

            th2 = rotate_right(socket, peg, newcm, config, pic)
            new_peg = Kinematics.get_peg_from_config(peg, config)

			if pic
	            Draw.draw_sth(socket, (0,0,0), 3.)
	            # Draw.draw_sth(new_peg, (1,1,0), 2.)
	            # Draw.draw_sth(peg, (1,0,1), 2.)
			end

            if th1 >= th2
                push!(th1_set, th1)
                push!(th2_set, th2)
            else
                push!(th1_set, th2)
                push!(th2_set, th1)
            end


            push!(th, abs(th1-th2))

            push!(id_set, v)
        end
    end

    sort!(th1_set)
    sort!(th2_set)
    max_th = th1_set[length(th1_set)]-th2_set[1]
    return th1_set, th2_set, max_th

end

function capture(th::Float64, dx::Float64, peg::myType.myPeg, socket::mySocket, ep=0.1)
    n_p = length(peg.points)
	n_s = length(socket.points)
    result1 = false
    result2 = false

	ls = vector_length([socket.points[1].x - socket.points[n_s].x
						socket.points[1].y - socket.points[n_s].y])/2

	lp = vector_length([peg.points[1].x - peg.points[n_p].x
						peg.points[1].y - peg.points[n_p].y])/2


	if ls >= lp*cos(th) + dx - ep
		result1 = true
	end

	v1 = [socket.points[1].x - socket.points[2].x, socket.points[1].y - socket.points[2].y, 0]
	v2 = [socket.points[n_s].x - socket.points[n_s-1].x, socket.points[n_s].y - socket.points[n_s-1].y, 0]

	v3 = [peg.points[1].x - socket.points[2].x, peg.points[1].y - peg.points[2].y, 0]
	v4 = [peg.points[n_p].x - socket.points[n_p-1].x, peg.points[n_p].y - peg.points[n_p-1].y, 0]




	th_op = compute_abs_angle(v1, v2)/2
	th_peg = compute_abs_angle(v3, v4)/2

	if th_op >= th_peg+th
		result2 = true
	end



	return result1||result2
end

function capture_false(th::Float64, dx::Float64, peg::myType.myPeg, socket::mySocket, ep=0.1)
    n = length(peg.points)

    result1 = false
    result2 = false

    if iseven(n)
		mid = Int(n/2)
        v1 = [peg.points[mid-1].x - peg.points[mid].x, peg.points[mid-1].y - peg.points[mid].y, 0]
        v2 = [socket.points[1].x - socket.points[2].x, socket.points[1].y - socket.points[2].y, 0]

        th1 = compute_signed_angle(v1, v2)

		mid = mid+1
        v3 = [peg.points[mid+1].x - peg.points[mid].x, peg.points[mid+1].y - peg.points[mid].y, 0]
        v4 = [socket.points[length(socket.points)].x - socket.points[length(socket.points)-1].x,
        socket.points[length(socket.points)].y - socket.points[length(socket.points)-1].y, 0]
        th2 = compute_signed_angle(v3, v4)

        if th1+th>=0 && th2-th<=0
            result1 = true
        end
    else
        mid = Int((n+1)/2)
        v1 = [peg.points[mid-1].x - peg.points[mid].x, peg.points[mid-1].y - peg.points[mid].y, 0]
        v2 = [socket.points[1].x - socket.points[2].x, socket.points[1].y - socket.points[2].y, 0]

        th1 = compute_signed_angle(v1, v2)

        v3 = [peg.points[mid+1].x - peg.points[mid].x, peg.points[mid+1].y - peg.points[mid].y, 0]
        v4 = [socket.points[length(socket.points)].x - socket.points[length(socket.points)-1].x,
        socket.points[length(socket.points)].y - socket.points[length(socket.points)-1].y, 0]
        th2 = compute_signed_angle(v3, v4)

        if th1+th>=0 && th2-th<=0
            result1 = true
        end
    end

    if socket.points[2].x >= -dx && socket.points[length(socket.points)-1].x <= dx
        result2 = true
    end

    return result1 && result2

end

function v_ratio(v1, v2)
	l = vector_length(v1)
	return abs((v1[1]*v2[1]+v1[2]*v2[2])/l^2)

end



function compute_t(peg::myType.myPeg, socket::myType.mySocket)
	n_s = length(socket.edges)
	n_p = length(peg.points)
	p_mid = socket.points[2]
	if n_s > 6 || n_s < 3 || n_p > 5 || n_p < 3
		error("no such joint")
	elseif n_s == 4 && n_p ==3
		p_i = socket.points[2]
		p_j = socket.points[1]
		p_p = peg.points[1]
		v1 = [p_j.x - p_i.x, p_j.y - p_i.y]
		v2 = [p_p.x - p_i.x, p_p.y - p_i.y]
		t = v_ratio(v1, v2)
		return t, t

	elseif n_s == 4 && n_p == 4
		p_i = socket.points[2]
		p_j = socket.points[1]
		p_p1 = peg.points[1]
		p_p2 = peg.points[2]
		v1 = [p_j.x - p_i.x, p_j.y - p_i.y]
		v2 = [p_p1.x - p_i.x, p_p.y - p_i.y]
		v3 = [p_j.x - p_i.x, p_j.y - p_i.y]
		v4 = [p_p2.x - p_i.x, p_p.y - p_i.y]
		t_up = v_ratio(v1, v2)
		t_low = v_ratio(v3, v4)
		return t_up, t_low

	elseif n_s > 4  && n_p > 3
		p_i = socket.points[2]
		p_j1 = socket.points[1]
		p_j2 = socket.points[3]
		p_p1 = peg.points[1]
		p_p2 = peg.points[2]
		v1 = [p_j1.x - p_i.x, p_j1.y - p_i.y]
		v2 = [p_p1.x - p_i.x, p_p1.y - p_i.y]
		v3 = [p_j2.x - p_i.x, p_j2.y - p_i.y]
		v4 = [p_p2.x - p_i.x, p_p2.y - p_i.y]
		t_up = v_ratio(v1, v2)
		t_low = v_ratio(v3, v4)
		return t_up, t_low
	end


end

function ep_trans_t(peg::myType.myPeg, socket::myType.mySocket, ep=0.1)
	p1 = socket.points[1]
	p2 = socket.points[2]
	v1 = [p2.x - p1.x, p2.y - p1.y, 0]
	th1 = compute_abs_angle(v1, [1,0,0])
	l1 = vector_length(v1)
	lep1 = ep/tan(th1)
	t1 = lep1/l1

	p3 = socket.points[3]
	v2 = [p2.x - p3.x, p2.y - p3.y, 0]
	th2 = compute_abs_angle(v2, [1,0,0])
	l2 = vector_length(v2)
	lep2 = ep/tan(th2)
	t2 = lep2/l2


	return t1, t2

end


function test_ad(i)
	peg, socket = Construct_Peg_Socket.read_from_file(i)
	peg = Adjust.get_adj_peg(peg, socket, 0.01)


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


function get_capture(peg, socket)
    th = pi/6
    dx = abs(socket.points[length(socket.points)].x - socket.points[1].x)/6
    return capture(th, dx, peg, socket)

end

function concave_dis(socket::myType.mySocket, ep = 0.1)
	n = length(socket.points)
	if n != 5 && n!=6
		error("wrong concave_dis")
	end
	return (distance(socket.points[2], socket.points[n-1]) - 2ep)
end

function peg_dis(peg::myType.myPeg)
	n = length(peg.points)
	if n != 4 && n!=5
		error("wrong peg_dis")
	end

	return distance(peg.points[2], peg.points[n-1])
end

function get_joint(peg::myType.myPeg,socket::myType.mySocket, F = [0,0,0],
	convex = false, svf = false)
	# To get joint, we first optimize peg with respect to t under constrains of capture and limit t
	# Then optimize the socket by rotating up and low edges of socket along the direction that reduce MRA
	# The low socket edge has constraints of d <= D
	# The up socket edge has constraints of caoptrue and not breaking sink

	Insertion_closed = false
	Stability_closed = true


	#svf = true



	ori_peg = deepcopy(peg)
	ori_socket = deepcopy(socket)


	peg, socket, m = optimizing_peg_and_capture(peg, socket, F, convex, svf)


	thset, fth, sth= Rotation.get_rotate_th(peg, socket, F, Stability_closed, convex)

	th_mid = fth[1] + sth[1]

	fname = "Data.txt"
	open(fname,"a") do io

		 println(io,"th_mid =",th_mid)

	end

	sink = Insertion.get_sink(peg, socket, F, Insertion_closed, convex)

	n_sink = length(sink)

	dth = 0.005

	th_limit = Insertion_gradient.limit_th(ori_peg, ori_socket, F, Insertion_closed, convex)

	rotate_th = edge_rotate_th(ori_socket, socket)


	k_set =[]


	while length(sink) > 0
		# println("here7")
		g = Insertion.get_graph_edges(peg, socket, F, Insertion_closed, convex)

		g1 = g
		count = 1
		new_peg = deepcopy(peg)
		new_socket = deepcopy(socket)

		while count <=100 && g1 == g
			# println("here5")
			r, index, k = Stability_gradient.socket_gradient(peg, socket, F, convex, dth)


			push!(k_set, r, index, k)


			if length(k_set) >= 3
				deleteat!(k_set, 1)
			end

			if length(k_set) == 2

					if k_set[1][1] == k_set[2][1] && k_set[1][2] == k_set[2][2] && k_set[1][3] * k_set[2][3] < 0

						dth = dth*1.5

						r, index, k = Stability_gradient.socket_gradient(peg, socket, F, convex, dth)

						pop!(k_set)
						push!(k_set, r, index, k)

					end

			end

			if k == nothing || k == 0 || dth > 0.05
				return peg, socket
			end

			socket_new = Rotation.rotate_socket(socket, r, index, k*dth)
			peg_new = Adjust.get_adj_peg(peg, socket_new, 0.)

			g1 = Insertion.get_graph_edges(peg_new, socket_new, F, Insertion_closed, convex)

			if g1 == g
				socket = socket_new
				peg = peg_new
				new_peg = peg_new
				new_socket = socket_new
			else
				new_peg = peg_new
				new_socket = socket_new
				break
			end

			count += 1
			PyPlot.clf()
			PyPlot.axis("off")
			Draw.draw_sth(ori_peg, (0,0,1), 3.)
			Draw.draw_sth(ori_socket, (0,0,0), 3.)
			Draw.draw_sth(peg, (1,0,0), 3.)
			Draw.draw_sth(socket, (0.5,.4,.2), 3.)

			# println("count:", count)
			m +=1
			if svf
				filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
				PyPlot.savefig(filename)
			end
		end

		# println("here4")
		sink = Insertion.get_sink(new_peg, new_socket, F, Insertion_closed, convex)
		# println(sink)

		if length(sink) > 0
			peg = new_peg
			socket = new_socket
			count += 1
			# println("plot6")
			PyPlot.clf()
			PyPlot.axis("off")
			Draw.draw_sth(ori_peg, (0,0,1), 3.)
			Draw.draw_sth(ori_socket, (0,0,0), 3.)
			Draw.draw_sth(peg, (1,0,0), 3.)
			Draw.draw_sth(socket, (0.5,.4,.2), 3.)

			# println("count:", count)
			m +=1
			if svf
				filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
				PyPlot.savefig(filename)
			end

		else
			# println("sink change")
			 # peg = new_peg
			 # socket = new_socket
			 break
		end

	end









	return peg, socket




end

function edge_rotate_th(socket::myType.mySocket, new_socket::myType.mySocket)

	n = length(socket.points)
	v1 = [socket.points[2].x - socket.points[1].x,
	socket.points[2].y - socket.points[1].y,
	0]

	v2 = [new_socket.points[2].x - new_socket.points[1].x,
	new_socket.points[2].y - new_socket.points[1].y,
	0]

	return compute_signed_angle(v1, v2)
end

function test_k_peg(peg::myType.myPeg, socket::myType.mySocket,
	F, Stability_closed, convex, limit_idex, dt, id, sign)

	dt[id] = dt[id]*2



	k_peg = Stability_gradient.peg_gradient(peg, socket,
		F, Stability_closed, convex, limit_idex)

		# println("k_peg, id:", k_peg, ",", id)
		# println("sign:", sign)
	if k_peg[id] == -sign || dt[id] >= 0.5
		# println("return k_peg:", k_peg)
		return k_peg
	else

		k_peg = test_k_peg(peg, socket,
			F, Stability_closed, convex, limit_idex, dt, id, sign)
	end

end

function optimizing_peg_and_capture(peg::myType.myPeg,socket::myType.mySocket, F = [0,0,0],
	 convex = false, svf = false, d_t = 0.01)
	#optimize peg with respect to t under constrains of capture and limit t
	Insertion_closed = false
	Stability_closed = true



	n_t = Int(floor(length(peg.points)/2))
	d_t = 0.01
	dth = 0.005

	#println(n_t)
	dt = []
	for i = 1:n_t
		push!(dt, d_t)
	end


	ori_peg = deepcopy(peg)
	ori_socket = deepcopy(socket)

	sink = Insertion.get_sink(peg, socket, F, Insertion_closed, convex)

	th_limit = Insertion_gradient.limit_th(peg, socket, F, Insertion_closed, convex)
	rotate_th = 0

	n_sink = length(sink)
	m = 0

	limit_idex = []

	k_set = []
	test = 0

	if !Insertion.is_convex(socket)

		while length(sink) > 0
			# println("here1")
			new_peg = deepcopy(peg)
			new_socket = deepcopy(socket)

			g = Insertion.get_graph_edges(peg, socket, F, Insertion_closed, convex)
			count = 1
			g1 = g
			g2 = g


			while g1 == g && g2 == g && count <=100

				while get_capture(peg, socket) && count <=100

					full = Adjust.check_peg_limit(peg, socket)



					k_peg = Stability_gradient.peg_gradient(peg, socket,
						F, Stability_closed, convex, limit_idex)

					# push!(k_set, k_peg)
					#
					#
					# if length(k_set) >= 3
					# 	deleteat!(k_set, 1)
					# 	test = k_set[1] .* k_set[2]
					# end
					#
					#
					# for j =  1 : length(test)
					# 	if test[j] < 0
					# 		# peg_test1  = Adjust.get_adj_peg(peg, socket, k_peg.*dt)
					# 		# return peg, peg_test1, socket
					#
					# 		sign = k_peg[j]
					# 		kpeg = test_k_peg(peg, socket,
					# 			F, Stability_closed, convex, limit_idex, dt, j, sign)
					#
					# 		# pop!(k_set)
					# 		# push!(k_set, k_peg)
					# 		k_set = []
					# 		test = 0
					#
					# 	end
					# end

					# println(full, k_peg)
					# println(dt)

					if full || unique(k_peg) == [0]
						return peg, socket, m
					end


					peg1  = Adjust.get_adj_peg(peg, socket, k_peg.*dt)

					g1 = Insertion.get_graph_edges(peg1, socket, F, Insertion_closed, convex)


					if g1 == g && get_capture(peg1, socket)
						new_socket = socket
						peg = peg1
						new_peg = peg1
					elseif g1 != g
						new_peg = peg1
						new_socket = socket
						break
					elseif !get_capture(peg1, socket)
						peg = peg1
						new_peg = peg1
						new_socket = socket
						break
					end

					# println("plot1")

					PyPlot.clf()
					PyPlot.axis("off")
					Draw.draw_sth(ori_peg, (0,0,1), 3.)
					Draw.draw_sth(ori_socket, (0,0,0), 3.)
					Draw.draw_sth(peg, (1,0,0), 3.)
					Draw.draw_sth(socket, (0.5,.4,.2), 3.)
					PyPlot.pause(0.1)


					count += 1
					# println("count:", count)
					m += 1
					if svf
						filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
						PyPlot.savefig(filename)
					end
				end

				while !get_capture(peg, socket) && rotate_th < th_limit && count <=100
					socket2 = Rotation.capture_rotate(socket, dth)
					peg2 = Adjust.get_adj_peg(peg, socket, 0.)

					rotate_th = edge_rotate_th(ori_socket, socket2)


					if rotate_th >= th_limit
						return peg, socket, m
					end

					g2 = Insertion.get_graph_edges(peg2, socket2, F, Insertion_closed, convex)
					if g2 == g && !get_capture(peg2, socket2)
						socket = socket2
						peg = peg2
						new_socket = socket2
						new_peg = peg2
					elseif g2 != g
						new_socket = socket2
						new_peg = peg2
						break
					elseif get_capture(peg2, socket2)
						socket = socket2
						peg = peg2
						new_socket = socket2
						new_peg = peg2
						break
					end

					count += 1
					# println("plot2")
					PyPlot.clf()
					PyPlot.axis("off")
					Draw.draw_sth(ori_peg, (0,0,1), 3.)
					Draw.draw_sth(ori_socket, (0,0,0), 3.)
					Draw.draw_sth(peg, (1,0,0), 3.)
					Draw.draw_sth(socket, (0.5,.4,.2), 3.)
					PyPlot.pause(0.1)
					# println("count:", count)

					m+=1
					if svf
						filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
			        	PyPlot.savefig(filename)
					end
				end

				if g1 == g && g2 == g
					count += 1
					# println("plot3")
					PyPlot.clf()
					PyPlot.axis("off")
					Draw.draw_sth(ori_peg, (0,0,1), 3.)
					Draw.draw_sth(ori_socket, (0,0,0), 3.)
					Draw.draw_sth(peg, (1,0,0), 3.)
					Draw.draw_sth(socket, (0.5,.4,.2), 3.)
					PyPlot.pause(0.1)
					# println("count:", count)

					m+=1
					if svf
					filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
					PyPlot.savefig(filename)
					end
					continue
				else
					break
				end
			end



			# println("here4")
			sink = Insertion.get_sink(new_peg, new_socket, F, Insertion_closed, convex)
			# println(sink)

			if length(sink) > 0
				peg = new_peg
				socket = new_socket

				# println("plot4")

				PyPlot.clf()
				PyPlot.axis("off")
				Draw.draw_sth(ori_peg, (0,0,1), 3.)
				Draw.draw_sth(ori_socket, (0,0,0), 3.)
				Draw.draw_sth(peg, (1,0,0), 3.)
				Draw.draw_sth(socket, (0.5,.4,.2), 3.)
				PyPlot.pause(0.1)

				# println("count:", count)
				m += 1
				if svf
					filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
					PyPlot.savefig(filename)
				end
			else

				# println("sink change")
				return new_peg, new_socket, m
				 # peg = new_peg
				 # socket = new_socket
				 break
			end

		end

		return peg, socket, m

	else

		while length(sink) > 0
			# println("here1")
			new_peg = deepcopy(peg)
			new_socket = deepcopy(socket)

			g = Insertion.get_graph_edges(peg, socket, F, Insertion_closed, convex)
			count = 1
			g1 = g
			g2 = g


			while g1 == g && g2 == g && count <=100

				while get_capture(peg, socket) && count <=100
					# println("here 6")

					full = Adjust.check_peg_limit(peg, socket)



					k_peg = Stability_gradient.peg_gradient(peg, socket,
						F, Stability_closed, convex, limit_idex)

					peg2  = Adjust.get_adj_peg(peg, socket, k_peg.*dt)

					if !get_capture(peg2, socket)
						push!(limit_idex , 1)
					end

					k_peg = Stability_gradient.peg_gradient(peg, socket,
						F, Stability_closed, convex, limit_idex)



					push!(k_set, k_peg)

					if length(k_set) >= 3
						deleteat!(k_set, 1)
						test = k_set[1] .* k_set[2]
					end


					for j =  1 : length(test)
						if test[j] < 0
							# peg_test1  = Adjust.get_adj_peg(peg, socket, k_peg.*dt)
							# return peg, peg_test1, socket
							dt[j] -= d_t/10
							if dt[j] <= d_t/4
								push!(limit_idex , j)
							end
						end
					end

					# println(full, k_peg)
					# println(dt)

					if full || unique(k_peg) == [0]
						return peg, socket, m
					end


					peg1  = Adjust.get_adj_peg(peg, socket, k_peg.*dt)


					g1 = Insertion.get_graph_edges(peg1, socket, F, Insertion_closed, convex)


					if g1 == g && get_capture(peg1, socket)
						peg = peg1
						new_peg = peg1
					elseif g1 != g
						new_peg = peg1
						break
					elseif !get_capture(peg1, socket)
						return peg, socket, m
					end

					# println("plot1")

					PyPlot.clf()
					PyPlot.axis("off")
					Draw.draw_sth(ori_peg, (0,0,1), 3.)
					Draw.draw_sth(ori_socket, (0,0,0), 3.)
					Draw.draw_sth(peg, (1,0,0), 3.)
					Draw.draw_sth(socket, (0.5,.4,.2), 3.)
					PyPlot.pause(0.1)


					count += 1
					# println("count:", count)
					m += 1
					if svf
						filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
						PyPlot.savefig(filename)
					end
				end

				if !get_capture(peg, socket)
					# println("cannot be captured")

					count += 1
					# println("plot5")
					PyPlot.clf()
					PyPlot.axis("off")
					Draw.draw_sth(ori_peg, (0,0,1), 3.)
					Draw.draw_sth(ori_socket, (0,0,0), 3.)
					Draw.draw_sth(peg, (1,0,0), 3.)
					Draw.draw_sth(socket, (0.5,.4,.2), 3.)
					PyPlot.pause(0.1)
					# println("count:", count)

					m+=1
					if svf
					filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
					PyPlot.savefig(filename)
					end
					return peg, socket, m
				end

				if g1 == g && g2 == g
					count += 1
					# println("plot3")
					PyPlot.clf()
					PyPlot.axis("off")
					Draw.draw_sth(ori_peg, (0,0,1), 3.)
					Draw.draw_sth(ori_socket, (0,0,0), 3.)
					Draw.draw_sth(peg, (1,0,0), 3.)
					Draw.draw_sth(socket, (0.5,.4,.2), 3.)
					PyPlot.pause(0.1)
					# println("count:", count)

					m+=1
					if svf
					filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
					PyPlot.savefig(filename)
					end
					continue
				else
					break
				end



			end



			# println("here4")
			sink = Insertion.get_sink(new_peg, new_socket, F, Insertion_closed, convex)
			# println(sink)

			if length(sink) > 0
				peg = new_peg
				socket = new_socket

				# println("plot4")

				PyPlot.clf()
				PyPlot.axis("off")
				Draw.draw_sth(ori_peg, (0,0,1), 3.)
				Draw.draw_sth(ori_socket, (0,0,0), 3.)
				Draw.draw_sth(peg, (1,0,0), 3.)
				Draw.draw_sth(socket, (0.5,.4,.2), 3.)
				PyPlot.pause(0.1)

				# println("count:", count)
				m += 1
				if svf
					filename = "./anime/" * string(m, base = 10, pad = 4) * ".png"
					PyPlot.savefig(filename)
				end
			else

				# println("sink change")
				return new_peg, new_socket, m
				 # peg = new_peg
				 # socket = new_socket
				 break
			end

		end

		return peg, socket, m
	end


end


function check_peg(i)

	closed = true
	F = [0,0,0]
    peg1, socket = Construct_Peg_Socket.read_from_file(i)
	peg = Adjust.get_adj_peg(peg1, socket, 0.)
	convex = Insertion.is_convex(socket)
	#get_rotate_th(peg,socket, F, closed, convex)
	peg1, peg2, socket = optimizing_peg_and_capture(peg,socket)

	thset1, fth1, sth1 = Rotation.get_rotate_th(peg1, socket, F, closed, convex, true)
	thset2, fth2, sth2 = Rotation.get_rotate_th(peg2, socket, F, closed, convex, true)
	# println("maxth1: ", [fth1, sth1])
	# println("maxth2: ", [fth2, sth2])
	# println("thset1:", thset1)
	# println("thset2:", thset2)



	# get_rotate_th(peg, socket, true)

end

function test_peg(i)
	Insertion_closed = false
	Stability_closed = true


	F = [0,0,0]
    peg1, socket = Construct_Peg_Socket.read_from_file(i)
	peg = Adjust.get_adj_peg(peg1, socket, 0.)
	convex = Insertion.is_convex(socket)
	#get_rotate_th(peg,socket, F, closed, convex)
	svf = false
	peg, socket, m = optimizing_peg_and_capture(peg, socket, F, convex, svf)

	#Insertion.model_graph(peg, socket,F, Insertion_closed, convex)

	# get_rotate_th(peg, socket, true)

end



function test(i)

	Insertion_closed = false
	Stability_closed = true


	F = [0,0,0]
    peg1, socket = Construct_Peg_Socket.read_from_file(i)
	peg = Adjust.get_adj_peg(peg1, socket, 0.)
	convex = Insertion.is_convex(socket)
	thset, fth, sth = Rotation.get_rotate_th(peg, socket, F, Stability_closed, convex)

	num_p = length(peg.points)
	num_e = length(socket.edges) - 1
	th_ori = fth[1] + sth[1]
	fname = "Data.txt"
	open(fname,"a") do io
         # println(io,"n, m =",[num_p, num_e])
		 # println(io,"th_ori =",th_ori)

    end


	peg, socket = get_joint(peg, socket, F, convex, true)
	thset, fth, sth = Rotation.get_rotate_th(peg, socket, F, Stability_closed, convex)
	# println("maxth1: ", fth)
	# println("maxth2: ", sth)
	# println(thset)

	th_final = fth[1] + sth[1]

	open(fname,"a") do io

		 # println(io,"th_final =",th_final)

	end

	jointfile = "joint.txt"


	open(jointfile,"a") do io

		 println(io,"n, m =", [num_p, num_e])

	end

	for  i = 1: length(peg.points)
		open(jointfile,"a") do io

			 println(io,"peg:[index, x, y] =", [i, peg.points[i].x, peg.points[i].y])

		end
	end

	for  j = 1: length(socket.points)
		open(jointfile,"a") do io

			 println(io,"socket:[index, x, y] =", [j, socket.points[j].x, socket.points[j].y])

		end
	end


	#Insertion.model_graph(peg, socket,F, Insertion_closed, convex)
	# get_rotate_th(peg, socket, true)

end

function test_adjsocket(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
	Draw.draw_sth(socket)
	socket = adj_upsocket(socket, 0.02)
	Draw.draw_sth(socket)
	socket = adj_lowsocket(socket, 0.02)
	Draw.draw_sth(socket)
end

function test_rotate(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
	Rotation.get_rotate_th(peg, socket, true)

end

function test_adpeg(i)
	peg,socket = Construct_Peg_Socket.read_from_file(i)
	peg = Adjust.get_adj_peg(peg, socket, 0)
	Draw.draw_sth(socket)
	Draw.draw_sth(peg)
end






end
