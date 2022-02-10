module Construct_Peg_Socket

import ..myType
import ..ConcaveTest

#construct the socket by given points and then construct the peg by inputed contacts
function read_from_file(i::Int)
	peg_file = "./joints_defined_by_location/peg-" * string(i) * ".dat"
	socket_file = "./joints_defined_by_location/socket-" * string(i) * ".dat"

	f_p = open(peg_file)
	lines = readlines(f_p)

	n = parse(Int, lines[1])

	peg = myType.myPeg()
	peg.points = Array{myType.Indexed_point}(undef, n)
	for i = 1:n
		result = split(lines[i+1], " ")
		peg.points[i] = myType.Indexed_point(i,
			parse(Float64, result[1]), parse(Float64, result[2]))

	end



	close(f_p)

	f_s = open(socket_file)
	lines = readlines(f_s)
	n = parse(Int, lines[1])

	socket = myType.mySocket()
	socket.points = Array{myType.Indexed_point}(undef, n)
	for i = 1:n
		result = split(lines[i+1], " ")
		socket.points[i] = myType.Indexed_point(i,
			parse(Float64, result[1]), parse(Float64, result[2]))
	end

	socket.edges = Array{myType.reference_edge}(undef, n-1)
	for i = 1:n-1
		next = (i+1)%n
		if next == 0
			next = n
		end

		socket.edges[i] = myType.reference_edge(i, i, next);
	end


	close(f_s)


	return peg, socket
end

function get_contact_ids(i)
	contacts_file = "./joints_defined_by_location/contacts-" * string(i) * ".dat"
	c_ids = Array{myType.contact_id}(undef, 0)
	f_c = open(contacts_file)
	c_lines = readlines(f_c)

	m = length(c_lines)
	for i = 1:m
		result = split(c_lines[i], " ")

		push!(c_ids, myType.contact_id(parse(Int, result[1]), parse(Int, result[2])))
	end
	close(f_c)
	return c_ids
end

function get_socket(points::Array{myType.Indexed_point})
	n = length(points)
	socket = myType.mySocket()
	socket.points = Array{myType.Indexed_point}(undef, n)
	socket.edges = Array{myType.reference_edge}(undef,n-1)
	for i = 1:n
		next = (i+1)%n
		if next == 0
			next = n
		end

		socket.points[i] = points[i]

		if i == n
			break
		end

		socket.edges[i] = myType.reference_edge(i, i, next);
	end
	return socket
end

function get_peg(socket::myType.mySocket, contacts::Array{myType.contact})
	n = length(contacts)
	peg = myType.myPeg()
	if !ConcaveTest.is_concave(socket)
		peg.points = Array{myType.Indexed_point}(undef, n)
		for i = 1:n
			e_id = contacts[i].e_id
			t = contacts[i].t
			p1 = socket.points[e_id]
			p2 = socket.points[e_id+1]
			location = get_contact_location(p1, p2, t)
			peg.points[i] = myType.Indexed_point(i, location[1], location[2])
		end
		return peg
	else
		sp_id = ConcaveTest.concave_points(socket)
		pp_id = deepcopy(sp_id)
		m = length(sp_id)
		for i = 1:m
			n1 = i-1
			n2 = 0
			for ct in contacts

				if ct.e_id < sp_id[i]
					n2 += 1
				end

			end

			pp_id[i] = n1+n2+1

		end

		peg.points = Array{myType.Indexed_point}(undef, n+m)

		count1 = 0

		for i =1:m+n
			if i in pp_id
				for j = 1:length(pp_id)
					if i == pp_id[j]
						peg.points[i] = socket.points[sp_id[j]]
						count1 += 1
					end
				end
			else
				count2 = i - count1

				e_id = contacts[count2].e_id
				t = contacts[count2].t
				p1 = socket.points[e_id]
				p2 = socket.points[e_id+1]
				location = get_contact_location(p1, p2, t)

				peg.points[i] = myType.Indexed_point(i, location[1], location[2])
			end

		end

	end
	return peg

end

function get_contact_location(p1::myType.Indexed_point, p2::myType.Indexed_point, t::Float64)
	v = [p2.x - p1.x, p2.y-p1.y]
	return [p1.x, p1.y] + t*v
end

function get_peg_and_socket(i)
	socket_file = "./joints_defined_by_contacts/socket-" * string(i) * ".dat"
	contacts_file = "./joints_defined_by_contacts/contacts-" * string(i) * ".dat"
	f_s = open(socket_file)
	s_lines = readlines(f_s)
	n = length(s_lines)
	ps = Array{myType.Indexed_point}(undef, 0)
	for i = 1:n
		result = split(s_lines[i], " ")
		push!(ps, myType.Indexed_point(i, parse(Float64, result[1]), parse(Float64, result[2])))
	end

	socket = get_socket(ps)

	close(f_s)

	f_c = open(contacts_file)
	c_lines = readlines(f_c)

	cs = Array{myType.contact}(undef, 0)
	m = length(c_lines)
	for i = 1:m
		result = split(c_lines[i], " ")

		push!(cs, myType.contact(parse(Int, result[1]), parse(Float64, result[2])))
	end

	peg = get_peg(socket, cs)

	close(f_c)

	return peg, socket
end

function get_contacts(i)
	contacts_file = "./joints_defined_by_contacts/contacts-" * string(i) * ".dat"
	f_c = open(contacts_file)
	c_lines = readlines(f_c)

	cs = Array{myType.contact}(undef, 0)
	m = length(c_lines)
	for i = 1:m
		result = split(c_lines[i], " ")

		push!(cs, myType.contact(parse(Int, result[1]), parse(Float64, result[2])))
	end
	close(f_c)
	return cs
end

function test(i)
	peg, socket = get_peg_and_socket(i)
	# println(peg)
	# println(socket)
end

end
