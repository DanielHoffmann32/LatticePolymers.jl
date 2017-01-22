module LatticePolymers

using StatPlots, StatsBase, DH32ParallelUtils

# package code goes here
export average_monomer_distance,
estimate_Boltzmann_weights,
initial_particles_placement,
MC_particles_around_polymer_1,
n_neighboring_particles,
n_procs_MC_particles_around_polymer_1,
process_MC_particles_around_polymer_1,
self_avoiding_cubic_lattice_random_walk,
self_avoiding_random_walk_in_box,
self_digest_without_attraction

"""
Input:

- number n of monomers (at least 3)

Output:

- 3D self-avoiding random walk on a cubic lattice with coordinates of monomers in an n x 3 array of integers [x_i, y_i, z_i]

"""
function self_avoiding_cubic_lattice_random_walk(n::Int64)
    if n<3
        error("n should be at least 3.")
    end
    Delta = [[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]]
    r = zeros(Int64, n, 3)
    r[2,:] = [1,0,0]
    i = 3
    while i<=n
        #add a new monomer in one of 5 random directions (away from previous monomer)
        delta = copy(Delta)
        dr = r[i-1,:] - r[i-2,:]
        for j in 1:6
            if delta[j]==-dr
                deleteat!(delta,j)
                break
            end
        end
        rnew = r[i-1,:] + delta[rand(1:5)]
            
        #check for overlap
        overlap = false
        for j in 1:(i-1)
            if rnew == r[j,:] #overlap -> remove monomer
                i = i-1
                overlap = true
                break
            end
        end
        if overlap == false #add monomer
            r[i,:] = rnew
            i = i+1
        end
    end
    r
end

"""
Input:

- n x 3 array r of integer coordinates of lattice polymer of length n

Ouput:

- average monomer distance in lattice units

"""
function average_monomer_distance(r::Array{Int64,2})
    n = length(r[:,1])
    s = 0.
    for i in 2:n
        for j in 1:(i-1)
            d = r[i,:]-r[j,:]
            s += sqrt(dot(d,d))
        end
    end
    2.*s/n/(n-1)
end

"""
We randomly place N_part particles on a cubic lattice box of length L. The state of each point in the box is marked by a value (elements of 3D float array "box"). If the value of a position in the box is marked with a value of zero or less it is free and can be used for the placement. Each successful placement at a position [x,y,z] is marked by adding a value m_part to that position of "box". 

Input:

- box: L^3 array of floats marking the states of all positions in the lattice box

- L: integer length of cubic lattice box

- N_part: number of particles to be placed

- m_part: number to be added at position x,y,z after successful placement there

Output:

- modified box

- positions of all placed particles (N_part x 3 integer array)
    
"""
function initial_particles_placement(box::Array{Float64,3}, L::Int64, N_part::Int64, m_part::Float64)
    r_part = zeros(Int64,N_part,3)

    if N_part > L^3
        error("more particles than volume elements")
    end

    #initial placement
    for part in 1:N_part
        placed = false
        while placed == false
            x,y,z=rand(1:L,3)
            if box[x,y,z]<=0. #every empty position (=0) or position with E<0 is considered
                box[x,y,z]+=m_part
                r_part[part,:] = [x,y,z]
                placed=true
            end
        end
    end
    
    box, r_part
end

"""
Computes the number of particles on a cubic lattice that are neighbors of a given particle "part", i.e. located on one of the six neighboring places (or less if part is located at the lattice boundary).

Input:

- part: index of given particle

- r_parts: n_parts x 3 array of integer particle positions on 3D lattice box

- box: lattice box (L x 3 float array) with all lattice position states marked

- mark: if a lattice position state is greater than "mark", this position counts as occupied

- L: length of lattice box

Output:

- number of neighbors of "part", i.e. number of occupied neighboring positions
    
"""
function n_neighboring_particles(part::Int64, r_parts::Array{Int64,2}, box::Array{Float64,3}, mark::Float64, L::Int64)
    x,y,z = r_parts[part,1:3]
    n_neigh=0
    if (x>1) && (box[x-1,y,z]>=mark)
        n_neigh+=1
    end
    if (x<L) && (box[x+1,y,z]>=mark)
        n_neigh+=1
    end
    if (y>1) && (box[x,y-1,z]>=mark)
        n_neigh+=1
    end
    if (y<L) && (box[x,y+1,z]>=mark)
        n_neigh+=1
    end
    if (z>1) && (box[x,y,z-1]>=mark)
        n_neigh+=1
    end
    if (z<L) && (box[x,y,z+1]>=mark)
        n_neigh+=1
    end
    n_neigh      
end

"""
We simulate an experiment with a box of enzymes digesting each other ("self-digest"). All enzymes are subjected to stochastic digestion by spatial neighbors (the more neighbors, the higher the probability of digestion), and they move in the box at random.

Input:

- Nsteps: number of digestion and move sweeps (per step we try for each enzyme a digestion and a move).

- box: lattice box (L x 3 float array) of lattice position states.

- L: length of lattice box.

- r_enz: Nenz x 3 integer array of enzyme positions on lattice

- N_enz: initial number of enzymes in the box

- m_enz: marker value for enzymes in lattice box (if box[x,y,z] > m_enz, then position x,y,z is occupied by an enzyme)

- p_surv: probability of an enzyme surviving the presence of a neighboring enzyme per step

Output:

- n_enzymes: integer array of length Nsteps, number of enzymes present at each of the Nstep time steps
    
"""
function self_digest_without_attraction(
    Nsteps::Int64, box::Array{Float64,3}, L::Int64, r_enz::Array{Int64,2}, Nenz::Int64,
    m_enz::Float64, p_surv::Float64)
    
    n_enzymes = zeros(Int64,Nsteps)

    for step in 1:Nsteps

        #proteolysis check and digest
        enz = Nenz
        while enz>0
            #Check neighborhoods of each enzyme for enzymes.
            n_neigh = n_neighboring_particles(enz, r_enz, box, m_enz, L)

            #The more enzymes in the direct neigborhood, the higher the probability of digestion.
            if n_neigh>0 
                if rand()>p_surv^n_neigh #digest -> remove 1 enzyme
                    box[r_enz[enz,1],r_enz[enz,2],r_enz[enz,3]]=0. #remove enzyme marker from box
                    r_enz = vcat(r_enz[1:(enz-1),:],r_enz[(enz+1):end,:]) #remove enzyme from array
                    Nenz -= 1
                end
            end
            enz -= 1
        end

        #try new placements
        for enz in 1:Nenz
            xtry,ytry,ztry=rand(1:L,3)
            x,y,z=r_enz[enz,:]
            box[x,y,z] = 0. #the old position is also an option
            if box[xtry,ytry,ztry]==0. #trial volume free
                box[xtry,ytry,ztry]=m_enz #move enz. there
                r_enz[enz,:] = [xtry,ytry,ztry]
            else
                box[x,y,z]=m_enz #failed move
            end
        end

        n_enzymes[step] = Nenz
    end
    n_enzymes
end

"""
We place a self-avoiding random walk polymer of length nmonos in a cubic lattice box of length L. Each box position occupied by a monomer will be marked by a value of m_poly. The non-occupied lattice positions neighboring a monomer are assigned an energy value of E.

Input:

- nmonos: number of monomers to be placed in the box.

- L: length of the cubic lattic box to be generated.

- m_poly: marker value for the presence of a monomer at box position.

- E: value of interaction energy assigned to lattice positions neighboring a monomer.

Output: box, r

- box: cubic lattice box (3D float array) with positions occupied by monomers marked (m_poly, E)

- r: nmonos x 3 integer array of lattice positions of monomers
    
"""
function self_avoiding_random_walk_in_box(nmonos::Int64, L::Int64, m_poly::Float64, E::Float64)
    if nmonos>L
        error("stretched polymer does not fit in box")
    end
    
    box = zeros(Float64,L,L,L)

    #put polymer into box
    r = self_avoiding_cubic_lattice_random_walk(nmonos)
    Lhalf = L/2
    cx = convert(Int64,round(mean(r[:,1])-Lhalf))
    cy = convert(Int64,round(mean(r[:,2])-Lhalf))
    cz = convert(Int64,round(mean(r[:,3])-Lhalf))
    for i in 1:nmonos
        r[i,:] -= [cx,cy,cz]
        box[r[i,1],r[i,2],r[i,3]]=m_poly
    end
    
    #fill all empty elements neighboring monomers with attractive energy values
    for i in 1:nmonos
        x = r[i,1]
        y = r[i,2]
        z = r[i,3]
        if (x < L) && (box[x+1,y,z]!=m_poly) #do this in each of the 6 directions
            box[x+1,y,z]+=E
        end
        if (x > 1) && (box[x-1,y,z]!=m_poly)
            box[x-1,y,z]+=E
        end
        if (y < L) && (box[x,y+1,z]!=m_poly)
            box[x,y+1,z]+=E
        end
        if (y > 1) && (box[x,y-1,z]!=m_poly)
            box[x,y-1,z]+=E
        end
        if (z < L) && (box[x,y,z+1]!=m_poly)
            box[x,y,z+1]+=E
        end
        if (z > 1) && (box[x,y,z-1]!=m_poly)
            box[x,y,z-1]+=E
        end
    end
    
    box, r
end

"""
MC simulation of N_steps sweeps of N_parts particles moving around a fixed polymer. If a monomer of the polymer and and a particle are lattice neighbors, they interact with an energy of E_part_poly_contact. The simulation is a Metropolis MC simulation at a value of RT of gas constant times temperature (same units as E_part_poly_contact). Particles are not allowed to overlap with other particles or polymer.

Input:

- N_steps: number of MC sweeps.

- box: cubic lattice box (3D float array) describing the state of each lattice position.

- L: length of lattice box.

- r_parts: N_parts x 3 integer array of positions of particles on lattice.

- N_parts: number of particles.

- m_part: marker value for a particle at a lattice box position.

- RT: gas constant times temperature in the same units as the interaction energy E_part_poly_contact.

- E_part_poly_contact: interaction energy between a monomer of the polymer and a neighboring particle.

Output: energies, particle_contacts

- energies: float array with total energy of system at each MC step

- particle_contacts: integer array with number of contacts between particles at each MC step.
    
"""
function MC_particles_around_polymer_1(
    N_steps::Int64, box::Array{Float64,3}, L::Int64, r_parts::Array{Int64,2}, N_parts::Int64, 
    m_poly::Float64, m_part::Float64, RT::Float64, E_part_poly_contact::Float64)

    if m_part <= m_poly
        error("m_part must be greater than m_poly")
    end 
    energies = zeros(N_steps)
    particle_contacts = zeros(Int64,N_steps)
    mark = m_part+7.*E_part_poly_contact #safe indicator for occupied element
    
    for step in 1:N_steps

        #try a new position for each particle
        for part in 1:N_parts
            xtry,ytry,ztry=rand(1:L,3)
            E_new = box[xtry,ytry,ztry]
            if E_new<=0. #space available => try move
                x,y,z=r_parts[part,:]
                E_old = box[x,y,z]-m_part
                if rand()<exp(-(E_new-E_old)/RT)
                    box[xtry,ytry,ztry]+=m_part #move part. there
                    box[x,y,z]-=m_part #remove it from the old place
                    r_parts[part,:] = [xtry,ytry,ztry]
                end
            end
        end
        
        #count number of particle pairs (pairs of particles that touch)
        n_neigh = 0
        E_total = 0.
        for part in 1:N_parts
            n_neigh += n_neighboring_particles(part, r_parts, box, mark, L)
            E_total += box[r_parts[part,1],r_parts[part,2],r_parts[part,3]]-m_part
        end
        particle_contacts[step] = n_neigh / 2
        energies[step] = E_total        
    end
    energies, particle_contacts
end

"""
A complete MC process run with multiple polymers, and for each polymer a MC run with particles. This program can be used to distribute processes over processors.

Input:

- N_polys: number of polymers (i.e. also MC runs) to be generated

- n_monos: number of monomers per polymer

- L: length of cubic lattice box

- burn_in: number of MC steps at the beginning to be discarded
    
- rest of input same as for MC_particles_around_polymer_1

Output: mean_monomer_dist, mean_energy, mean_contacts

- mean_monomer_dist: for each polymer the mean monomer distance

- mean_energy: for each MC run, sum of negative values in box (= interaction energies between polymer and particles)

- mean_contacts: mean number contacts between particle per MC run
"""
function process_MC_particles_around_polymer_1(
    N_polys::Int64, n_monos::Int64, L::Int64, N_steps::Int64, burn_in::Int64,
    m_poly::Float64, E_contact::Float64, N_parts::Int64, m_part::Float64, RT::Float64)

    mean_monomer_dist = zeros(N_polys)
    mean_energy = zeros(N_polys)
    mean_contacts = zeros(N_polys)

    for run in 1:N_polys
        box, r_poly = self_avoiding_random_walk_in_box(n_monos, L, m_poly, E_contact)
        box, r_parts = initial_particles_placement(box, L, N_parts, m_part)
        energies, enzyme_contacts = 
        MC_particles_around_polymer_1(
        N_steps, box, L, r_parts, N_parts, m_poly, m_part, RT, E_contact
        )
        mean_monomer_dist[run] = average_monomer_distance(r_poly)
        mean_energy[run] = mean(energies[burn_in:end])
        mean_contacts[run] = mean(enzyme_contacts[burn_in:end])
    end
    mean_monomer_dist, mean_energy, mean_contacts
end

"""
Distribute MC runs over processors in procIDs. Calls process_MC_particles_around_polymer_1 for each process and collects results.

Input:

- procIDs: int array of process IDs (produces with addprocs())

- rest of input as for process_MC_particles_around_polymer_1.

Output:

- mdist, mene, mcont as process_MC_particles_around_polymer_1

- additionally: array of Boltzmann weights of each element of the other arrays
    
"""
function n_procs_MC_particles_around_polymer_1(
                                               procIDs::Array{Int64,1},
                                               N_polys::Int64,
                                               n_monos::Int64,
                                               L::Int64,
                                               N_steps::Int64,
                                               burn_in::Int64,
                                               m_poly::Float64,
    E_contact::Float64,
    N_parts::Int64,
    m_part::Float64,
    RT::Float64)

    n_procs = length(procIDs)
    N_poly_i = split_int_lengths(N_polys, n_procs)
    refs = Array(Future, n_procs)
    for i_proc in 1:n_procs
        refs[i_proc] = @spawnat procIDs[i_proc] process_MC_particles_around_polymer_1(
                                                                                      N_poly_i[i_proc],
                                                                                      n_monos,
                                                                                      L,
                                                                                      N_steps,
                                                                                      burn_in,
        m_poly, E_contact, N_parts, m_part, RT)
    end

    mdist, mene, mcont = fetch(refs[1])
    for i_proc in 2:n_procs
        md,me,mc = fetch(refs[i_proc])
        append!(mdist,md)    
        append!(mene,me)    
        append!(mcont,mc)
    end
    
    mdist, mene, mcont, estimate_Boltzmann_weights(mene, RT)
end

"""
For a given sample of N energies (assumed to be representative), estimate the Boltzmann weights of each of the N states.

Input:

- energies: 1D array of energies (floats)

- RT: gas constant times temperature in same units as energy

Output:

- array of N Boltzmann weights

"""
function estimate_Boltzmann_weights(
                                     energies::Array{Float64,1},
                                     RT::Float64
                                     )
    N_energies = length(energies)
    Z = sum(exp(-energies[1:N_energies]/RT))
    exp(-energies[1:N_energies]/RT)/Z
end 

end # module
