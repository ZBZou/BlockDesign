module Visualization

import Base.Iterators

using PyCall

import PyPlot
import JSON
import LightGraphs: vertices, edges, src, dst
import MetaGraphs: props
import GeometryTypes
import GeometryTypes: Point2f0
import CoordinateTransformations
import LinearAlgebra

using Plots


import ..myType
import ..Curve
import ..VectorUtil
import ..buildMatrix


const MINT_GREEN = (144.0 / 255.0, 238.0 / 255.0, 144.0 / 255.0)
const WARM_GREY = (211.0 / 255.0, 211.0 / 255.0, 211.0 / 255.0)
const BLACK = (0.0, 0.0, 0.0)
const RED = (1.0, 0.0, 0.0)
const YELLOW = (1.0, 1.0, 0.0)
const WHITE = (1.0, 1.0, 1.0)
const BLUE = (0.0, 0.0, 1.0)

const matplotlib_patches = PyNULL()
const matplotlib_collections = PyNULL()




function __init__()
    PyPlot.matplotlib.use("Qt5Agg")

    copy!(matplotlib_patches, pyimport_conda("matplotlib.patches", "matplotlib"))
    copy!(matplotlib_collections, pyimport_conda("matplotlib.collections",
        "matplotlib"))
end


function convert_polygon_to_matplotlib_polygon(boundary_points)
    boundary_as_array = reduce(vcat, map(transpose, boundary_points))

    matplotlib_patches[:Polygon](boundary_as_array)
end

"""
Used to fix RGB color to RGBA array, if color is correct then return original
"""
function check_color(color_tuple)
    l_c = length(color_tuple)

    if !(l_c == 3 || l_c == 4)
        error("Incorrect color format")
    end

    return color_tuple
end


function draw_polygon(points::Array{myType.myPoint2d})
	# the points should be in form of a sequence of points
	# a closed shape, so, the first point should be appended
	# to the list, to draw
	p = points[1]
	push!(points, p)
	draw_line_seq(points)
end

function draw_line_seq(points::Array{myType.myPoint2d}, c=(1, 0, 0), lw=3.0)
	# a sequence of points, not necessarily a closed shape
	len = size(points, 1)
	# the number of points

	for i = 1:len - 1
		x1 = points[i].x
		y1 = points[i].y

		x2 = points[i+1].x
		y2 = points[i+1].y

		PyPlot.plot([x1, x2], [y1, y2], color=c, linewidth=lw)
	end

end

function convert_polygon_to_matplotlib_polygon(boundary_points)
    boundary_as_array = reduce(vcat, map(transpose, boundary_points))

    matplotlib_patches.Polygon(boundary_as_array)
end

function clear_axis_markers!(axis)
    PyPlot.box(on=nothing)
    # axis[:get_xaxis]()[:set_visible](false)
    # axis[:get_yaxis]()[:set_visible](false)
    axis.get_xaxis().set_visible(false)
    axis.get_yaxis().set_visible(false)
end


function draw_circle(patches, x::Float64, y::Float64, r::Float64;
    c=(1, 0, 0), f=true)
	# push!(patches, matplotlib_patches[:Circle]((x, y), r, color=c, fill=f))
	# println("f is: ", f)
	# if f
	# 	push!(patches, matplotlib_patches.Circle((x, y), r, color=c, fill=True))
	# else
	# 	push!(patches, matplotlib_patches.Circle((x, y), r, color=c, fill=False))
	# end

	push!(patches, matplotlib_patches.Circle((x, y), r, color=c, fill=f))

end

function setup_canvas()
	PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []

    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((Red..., alpha))


    # axis[:add_collection](collection)
    axis.add_collection(collection)

    return figure, axis
end

function draw_mat_polygon(patches, points)
	push!(patches, matplotlib_patches.Polygon(points))
end

function draw_patches()
	PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches1 = []

    draw_circle(patches, 1.0, 2.0, 1.0)
    draw_circle(patches1, 2.0, 3.0, 1.0)

    points = rand(5, 2)

    draw_mat_polygon(patches, points)


    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((RED..., alpha))
    collection.set_edgecolor(BLACK)
    collection1.set_color((WHITE..., 0))
    collection1.set_edgecolor(BLACK)


    # axis[:add_collection](collection)
    axis.add_collection(collection)
    axis.add_collection(collection1)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()
end


function draw_part(patches, curve::Curve.myCurve)
	points = Array{Float64}(undef, 0, 2)
	for i = 0:0.005:1
		x = Curve.curve_get_x(curve, i)
		points = vcat(points, LinearAlgebra.transpose(x))
	end

	# points are the points on the boundary of the curve
	# now, plow the points

	draw_mat_polygon(patches, points)

end

# input should be an array
# the contents in the array should be decided soon
function draw_fixels_as_circles(patches, contacts, curve::Curve.myCurve)
    for i = 1:size(contacts, 1)
        x, tangent, normal = Curve.curve_get(curve, contacts[i])

        cx = x[1] + normal[1] * 0.06
        cy = x[2] + normal[2] * 0.06

        draw_circle(patches, cx, cy, 0.04)
    end
end


function draw_wrenches(axis, patches, G)
	# A[:, 1]: get the first column of the matrix

	for i = 1:size(G, 2)
		wrench = G[:, i]

		# get the wrench

		if VectorUtil.vector_length(wrench, 2) < 0.001
			# draw a circle, pure torque
			draw_circle(patches, 0.0, 0.0, 0.05)
		else
			f = VectorUtil.planar_wrench_to_force(wrench)
			# draw a force at f[1], f[2]
			# vector: wrench[1]/2, wrench[2] / 2, 1
			# this use planar draw would be the best

			axis.arrow(f[1], f[2], wrench[1]/2, wrench[2]/2,
                head_width=0.1, head_length=0.1, fc="k", ec="k")
		end

	end

end


function test_draw_all()
	PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []

    curve = Curve.curve_ellipse_new(1.0, 1.0)
    # curve = Curve.curve_poly_load("./pawl.dat");

    draw_part(patches, curve)

    G = Matrix{Float64}(undef, 3, 3)
    G[1, 1] = 0.0
    G[1, 2] = 0.0
    G[1, 3] = 1.0
    G[2, 1] = -1.0
    G[2, 2] = 0.0
    G[2, 3] = 0.2
    G[3, 1] = 1.0
    G[3, 2] = 1.0
    G[3, 3] = 0.2

    draw_wrenches(axis, patches, G)

    contacts = [0.0, 0.3]
    mode = "rs"

    draw_fixels_as_circles(patches_fixels, contacts, curve)


    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((RED..., alpha))
    collection.set_edgecolor(BLACK)


    collection_fixels = matplotlib_collections.PatchCollection(patches_fixels)
    collection_fixels.set_color((BLUE..., alpha))
    collection_fixels.set_edgecolor(BLUE)

    # axis[:add_collection](collection)
    axis.add_collection(collection)
    axis.add_collection(collection_fixels)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()
end

function display_cone(curve::Curve.myCurve, contacts::Array{Float64}, G)

    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []

    draw_part(patches, curve)

    draw_wrenches(axis, patches, G)

    draw_fixels_as_circles(patches_fixels, contacts, curve)


    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((RED..., alpha))
    collection.set_edgecolor(BLACK)


    collection_fixels = matplotlib_collections.PatchCollection(patches_fixels)
    collection_fixels.set_color((BLUE..., alpha))
    collection_fixels.set_edgecolor(BLUE)

    # axis[:add_collection](collection)
    axis.add_collection(collection)
    axis.add_collection(collection_fixels)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()
end



function display_cone_union(curve::Curve.myCurve,
        contacts::Array{Float64}, cu::myType.ConeUnion)

    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []

    draw_part(patches, curve)

    # draw_wrenches(axis, patches, G)

    draw_fixels_as_circles(patches_fixels, contacts, curve)


    # loop to test cone

    xs = []
    ys = []
    total = 0
    stepsize = 0.001
    for s = 0:stepsize:1
        if buildMatrix.coneUnion_testpoint(cu, curve, s) == 0
            x, tangent, normal = Curve.curve_get(curve, s)
            if x[2] < 0
                if size(xs, 1) > 0
                    axis.plot(xs, ys, color="g")
                    xs = []
                    ys = []
                end
                continue
            end
            push!(xs, x[1]+normal[1]*0.02);
            push!(ys, x[2]+normal[2]*0.02);
            total += stepsize
        elseif size(xs, 1) > 0
            # println("here to draw")
            axis.plot(xs, ys, color="g")
            xs = []
            ys = []
        end
    end
    # plot
    println("total length: ", total)






    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    # collection.set_color((RED..., alpha))
    collection.set_color((MINT_GREEN..., alpha))
    collection.set_edgecolor(BLACK)


    collection_fixels = matplotlib_collections.PatchCollection(patches_fixels)
    collection_fixels.set_color((BLUE..., alpha))
    collection_fixels.set_edgecolor(BLUE)

    # axis[:add_collection](collection)
    axis.add_collection(collection)
    axis.add_collection(collection_fixels)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()


end


function display_polygon(points)
    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4

    draw_mat_polygon(patches, points);



    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((RED..., alpha))
    collection.set_edgecolor(BLACK)



    axis.add_collection(collection)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()

end

function add_fixel_socket(axis, width, depth)
    # currently, assuming a straight socket

    xs = []
    ys = []

    x_scale = width + 0.04
    y_scale = depth + 0.04

    push!(xs, -x_scale - 0.1)
    push!(xs, -x_scale + 0.2)
    push!(xs, x_scale - 0.2)
    push!(xs, x_scale + 0.1)


    push!(ys, 1)
    push!(ys, -y_scale)
    push!(ys, -y_scale)
    push!(ys, 1)


    axis.plot(xs, ys, linewidth=5, color="k", zorder=0);



end


function display_cone_with_socket(curve::Curve.myCurve,
        contacts::Array{Float64}, G, width, depth)

    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    add_fixel_socket(axis, width, depth)

    draw_part(patches, curve)

    draw_wrenches(axis, patches, G)

    draw_fixels_as_circles(patches_fixels, contacts, curve)


    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((MINT_GREEN..., alpha))
    collection.set_edgecolor(BLACK)


    collection_fixels = matplotlib_collections.PatchCollection(patches_fixels)
    collection_fixels.set_color((RED..., 1))
    collection_fixels.set_edgecolor(RED)

    # axis[:add_collection](collection)
    axis.add_collection(collection)
    axis.add_collection(collection_fixels)


    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()
end


# below are functions that may be used to generate figures for the paper





function plot_for_paper_polygon(curve::Curve.myCurve, width, depth)
    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    add_fixel_socket(axis, width, depth)

    draw_part(patches, curve)



    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((MINT_GREEN..., alpha))
    collection.set_edgecolor(BLACK)


    # collection_fixels = matplotlib_collections.PatchCollection(patches_fixels)
    # collection_fixels.set_color((RED..., 1))
    # collection_fixels.set_edgecolor(RED)

    # axis[:add_collection](collection)
    axis.add_collection(collection)
    # axis.add_collection(collection_fixels)


    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()

end



function add_fixel_socket_tilt(axis, width, depth)
    # currently, assuming a straight socket

    xs = []
    ys = []

    x_scale = width + 0.04
    y_scale = depth + 0.04

    push!(xs, -x_scale - 0.1 - 0.2)
    push!(xs, -x_scale + 0.2 - 0.32)
    push!(xs, x_scale - 0.2 - 0.32)
    push!(xs, x_scale + 0.1 - 0.22)


    push!(ys, 1 + 0.43)
    push!(ys, -y_scale + 0.43)
    push!(ys, -y_scale + 0.43)
    push!(ys, 1 + 0.43)


    axis.plot(xs, ys, linewidth=5, color="k", zorder=0);



end

function add_fixel_socket_multi(axis, width, depth)
    # currently, assuming a straight socket

    xs = []
    ys = []

    x_scale = width + 0.04
    y_scale = depth + 0.04

    # left top
    push!(xs, -x_scale - 0.35)
    # left mid
    push!(xs, -x_scale + 0.4)
    # left bottom
    push!(xs, -x_scale + 0.55)
    # right bottom
    push!(xs, x_scale - 0.55)
    # right mid
    push!(xs, x_scale - 0.4)
    # right top
    push!(xs, x_scale + 0.35)

    # left top
    push!(ys, 1)
    # left mid
    push!(ys, -y_scale/2)
    # left bottom
    push!(ys, -y_scale - 0.04)
    # right bottom
    push!(ys, -y_scale - 0.04)
    # right mid
    push!(ys, -y_scale/2)
    # right top
    push!(ys, 1)


    axis.plot(xs, ys, linewidth=5, color="k", zorder=0);



end

function plot_for_paper_tilt_polygon(curve::Curve.myCurve, width, depth)
    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    # add_fixel_socket_tilt(axis, width, depth)
    add_fixel_socket_multi(axis, width, depth)

    draw_part(patches, curve)



    # collection = matplotlib_collections[:PatchCollection](patches)
    collection = matplotlib_collections.PatchCollection(patches)
    # collection1 = matplotlib_collections.PatchCollection(patches1)
    # collection[:set_color]((RED..., alpha))
    collection.set_color((MINT_GREEN..., alpha))
    collection.set_edgecolor(BLACK)


    # collection_fixels = matplotlib_collections.PatchCollection(patches_fixels)
    # collection_fixels.set_color((RED..., 1))
    # collection_fixels.set_edgecolor(RED)

    # axis[:add_collection](collection)
    axis.add_collection(collection)
    # axis.add_collection(collection_fixels)


    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()

end

function draw_dialte_lines(points1, points2)


    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    # add_fixel_socket_tilt(axis, width, depth)
    xs1 = []
    ys1 = []
    xs2 = []
    ys2 = []

    for i = 1:size(points1, 1)
        push!(xs1, points1[i].x)
        push!(ys1, points1[i].y)
        push!(xs2, points2[i].x)
        push!(ys2, points2[i].y)
    end

    push!(xs1, points1[1].x)
    push!(ys1, points1[1].y)
    push!(xs2, points2[1].x)
    push!(ys2, points2[1].y)

    axis.plot(xs1, ys1, linewidth=5, color="k", zorder=0);
    axis.plot(xs2, ys2, linewidth=5, color="g", zorder=0);


    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()

end




function draw_dialte_lines_with_joint(points1, points2, x, y)


    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    # add_fixel_socket_tilt(axis, width, depth)
    xs1 = []
    ys1 = []
    xs2 = []
    ys2 = []

    for i = 1:size(points1, 1)
        push!(xs1, points1[i].x)
        push!(ys1, points1[i].y)
        push!(xs2, points2[i].x)
        push!(ys2, points2[i].y)
    end

    push!(xs1, points1[1].x)
    push!(ys1, points1[1].y)
    push!(xs2, points2[1].x)
    push!(ys2, points2[1].y)

    axis.plot(xs1, ys1, linewidth=5, color="k", zorder=0);
    axis.plot(xs2, ys2, linewidth=5, color="g", zorder=0);

    xs = []
    ys = []
    for i =1:size(x, 1)
        push!(xs, x[i])
        push!(ys, y[i])
    end

    push!(xs, x[1])
    push!(ys, y[1])

    axis.plot(xs, ys, linewidth=1, color="r", zorder=0)

    # println("xs: ", xs)
    # println("ys: ", ys)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()

end




function draw_dialte_lines_with_joint_and_rotate(points1, points2, x, y, nx, ny)


    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    # add_fixel_socket_tilt(axis, width, depth)
    xs1 = []
    ys1 = []
    xs2 = []
    ys2 = []

    for i = 1:size(points1, 1)
        push!(xs1, points1[i].x)
        push!(ys1, points1[i].y)
        push!(xs2, points2[i].x)
        push!(ys2, points2[i].y)
    end

    push!(xs1, points1[1].x)
    push!(ys1, points1[1].y)
    push!(xs2, points2[1].x)
    push!(ys2, points2[1].y)

    axis.plot(xs1, ys1, linewidth=5, color="k", zorder=0);
    axis.plot(xs2, ys2, linewidth=5, color="g", zorder=0);

    xs = []
    ys = []
    for i =1:size(x, 1)
        push!(xs, x[i])
        push!(ys, y[i])
    end

    push!(xs, x[1])
    push!(ys, y[1])

    push!(nx, nx[1])
    push!(ny, ny[1])

    axis.plot(xs, ys, linewidth=3, color="r", zorder=0)
    axis.plot(nx, ny, linewidth=3, color="b", zorder=0)

    # println("xs: ", xs)
    # println("ys: ", ys)

    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()

end


function display_intersection(points_error, points_rotate)
    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    # add_fixel_socket_tilt(axis, width, depth)
    xs1 = []
    ys1 = []
    xs2 = []
    ys2 = []

    for i = 1:size(points_error, 1)
        push!(xs1, points_error[i].x)
        push!(ys1, points_error[i].y)
        push!(xs2, points_rotate[i, 1])
        push!(ys2, points_rotate[i, 2])
    end

    push!(xs1, points_error[1].x)
    push!(ys1, points_error[1].y)
    push!(xs2, points_rotate[1, 1])
    push!(ys2, points_rotate[1, 2])

    axis.plot(xs1, ys1, linewidth=5, color="k", zorder=0);
    axis.plot(xs2, ys2, linewidth=5, color="g", zorder=0);


    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()
end


function display_peg_socket(socket::myType.mySocket, points::Array{myType.Indexed_point})
    PyPlot.close("all")
    figure, axis = PyPlot.subplots(figsize=(10, 10))
    # axis[:set_aspect]("equal")
    axis.set_aspect("equal")
    Visualization.clear_axis_markers!(axis)

    patches = []
    alpha = 0.4
    patches_fixels = []


    # add_fixel_socket_tilt(axis, width, depth)
    xs1 = []
    ys1 = []
    xs2 = []
    ys2 = []

    for i = 1:size(socket.points, 1)
        push!(xs1, socket.points[i].x)
        push!(ys1, socket.points[i].y)
    end


    for i = 1:size(points, 1)
        push!(xs2, points[i].x)
        push!(ys2, points[i].y)
    end

    push!(xs2, points[1].x)
    push!(ys2, points[1].y)


    axis.plot(xs1, ys1, linewidth=5, color="k", zorder=0);
    axis.plot(xs2, ys2, linewidth=5, color="g", zorder=0);


    axis.set_xlim((-5, 5))
    axis.set_ylim((-5, 5))

    # return figure, axis
    PyPlot.show()
end



function draw_sth(sth, c=(1, 0, 0), lw=2.0)


    p = Array{myType.Indexed_point}(undef, 0)

    for i = 1:length(sth.points)
        push!(p, sth.points[i])
    end

    push!(p, sth.points[1])



    #draw_polygon(p,c,lw)
    draw_line_seq(p, c, lw)

    PyPlot.show()

end

function draw_open_sth(sth, c=(1, 0, 0), lw=2.0)
    p = Array{myType.Indexed_point}(undef, 0)

    for i = 1:length(sth.points)
        push!(p, sth.points[i])
    end





    #draw_polygon(p,c,lw)
    draw_line_seq(p, c, lw)

    PyPlot.show()

end



function draw_line_seq(points::Array{myType.Indexed_point}, c=(1, 0, 0), lw=2.0)
	# a sequence of points, not necessarily a closed shape
	len = size(points, 1)
	# the number of points


	for i = 1:len - 1
		x1 = points[i].x
		y1 = points[i].y

		x2 = points[i+1].x
		y2 = points[i+1].y



		PyPlot.plot([x1, x2], [y1, y2], color=c, linewidth=lw)


	end



end

function draw_mlti_sth(socket::Array{myType.mySocket},peg, c=(1, 0, 0), lw=2.0)
    PyPlot.figure()
    for i = 1:length(socket)
        p = Array{myType.Indexed_point}(undef, 0)

        for j = 1:length(socket[i].points)
            push!(p, socket[i].points[j])
        end

        draw_line_seq(p, (1,1,0), lw)
    end

    p = Array{myType.Indexed_point}(undef, 0)

    for i = 1:length(peg.points)
        push!(p, peg.points[i])
    end

    push!(p, peg.points[1])



    #draw_polygon(p,c,lw)
    draw_line_seq(p, c, lw)

    PyPlot.show()
end

end # module
