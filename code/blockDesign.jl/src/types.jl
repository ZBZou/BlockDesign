module myType


struct myPoint2d
	x::Float64
	y::Float64
	myPoint2d() = new()
	myPoint2d(x, y) = new(x, y)
end

struct myPoint3d
	x::Float64
	y::Float64
	z::Float64
end



mutable struct Indexed_point
	index::Int
	x::Float64
	y::Float64

	function Indexed_point(point::Indexed_point)
		new(point.index, point.x, point.y)
	end

	function Indexed_point(id::Int, lx::Float64, ly::Float64)
		new(id, lx, ly)
	end

	function Indexed_point()
		new(0, 0, 0)
	end


end


mutable struct Indexed_edge
	index::Int
	s::Indexed_point
	e::Indexed_point

	function Indexed_edge(id, source::Indexed_point, target::Indexed_point)
		s = Indexed_point(source.index, source.x, source.y)
		e = Indexed_point(target.index, target.x, target.y)
		new(id, s, e)
	end

	function Indexed_edge()
		new(0, Indexed_point(), Indexed_point())
	end

end

mutable struct reference_edge
	index::Int
	s::Int
	e::Int
	reference_edge() = new()
	reference_edge(index, s, e) = new(index, s, e)
end

mutable struct contact_pair
	id::Int
	c_p::Indexed_point
	c_e::Indexed_edge

	# contact_pair(id, c_p, c_e) = new(id, c_p, c_e)
	function contact_pair(id::Int, point::Indexed_point, edge::Indexed_edge)
		# c_p = Indexed_point(point)
		# c_e = edge
		new(id, point, edge,)
	end

	function contact_pair()
		new(0, Indexed_point(), Indexed_edge())
	end
end

struct c_pair
	point::Int
	edge::Int
	mode::Int


	c_pair() = new()
	c_pair(point, edge, mode) = new(point, edge, mode)
end

mutable struct contact_mode
	cps::Array{contact_pair}

	function contact_mode(c_p::Array{contact_pair})
		new(c_p)
	end

end

mutable struct ConeUnion
  num_cones::Int
  F::Array{Matrix{Float64}}
end

mutable struct myPeg
	points::Array{Indexed_point}

	myPeg() = new()
	function myPeg(points::Array{Indexed_point})
		new(points)
	end

end

mutable struct mySocket
	points::Array{Indexed_point}
	edges::Array{reference_edge} # the references are to the index in points array
	mySocket() = new()
		function mySocket(points::Array{Indexed_point})

			n = length(points)
			edges = Array{reference_edge}(undef, n)
			for i = 1:n
				next = (i+1)%n
				if next == 0
					next = n
				end
				edges[i] = myType.reference_edge(points[i].index, points[i].index, points[next].index);
			end
			new(points, edges)
		end

end




mutable struct mySocketEdges
	edges::Array{Indexed_edge}

	function mySocketEdges(socket::mySocket)
		edges = Array{Indexed_edge}(undef, 0)
		for i = 1:length(socket.edges)
			push!(edges, Indexed_edge(i, socket.points[socket.edges[i].s], socket.points[socket.edges[i].e]))
		end
		new(edges)
	end

end

mutable struct myPegEdges
	edges::Array{Indexed_edge}
	function myPegEdges(peg::myPeg)
		edges = Array{Indexed_edge}(undef, 0)
		for i = 1:length(peg.edges)
			push!(edges, Indexed_edge(i, peg.points[peg.edges[i].s], peg.points[peg.edges[i].e]))
		end
		new(edges)
	end

end


mutable struct configSE2
	x::Float64
	y::Float64
	th::Float64

	configSE2() = new()
	configSE2(x, y, th) = new(x, y, th)

end

mutable struct contact
	e_id::Int
	t::Float64

	contact() = new()
	contact(e_id, t) = new(e_id, t)
end

mutable struct contact_id
	p_id::Int
	e_id::Int

	contact_id() = new()
	contact_id(p_id, e_id) = new(p_id, e_id)
end

mutable struct cm_id
	cps::Array{contact_id}

	cm_id() = new()
	cm_id(cps) = new(cps)

end


export myPoint2d, myPoint3d, Indexed_point, Indexed_edge, reference_edge
export contact_pair, contact_mode, c_pair
export ConeUnion
export myPeg, mySocket, mySocketEdges
export configSE2, contact, contact_id, cm_id





end # modulecm1::myType.contact_mode, cm2::myType.contact_mode
