using LatticePolymers
using Base.Test

procIDs = addprocs(2)

# write your own tests here

### average_monomer_distance ###
# a "polymer" of 4 monomers, with monomer i=1,2,3,4 at position (i,i,i) 
r = [1 1 1
     2 2 2
     3 3 3
     4 4 4]
@test isapprox(LatticePolymers.average_monomer_distance(r), 2.886751345948129)

### initial_particles_placement ###
L = 20
box = zeros(L,L,L)
N_part = 100
m_part = 1.
box, r = LatticePolymers.initial_particles_placement(box, L, N_part, m_part)
@test size(r) == (N_part, 3)
@test maximum(r) <= L
@test minimum(r) >= 1
@test isapprox(sum(box), N_part)
s = sum([box[r[i,1],r[i,2],r[i,3]] for i in 1:N_part])
@test isapprox(s, N_part)

### n_neighboring_particles ###
L = 12
box = zeros(L,L,L)
mark = 1.
r_parts = [1 1 1
           1 1 10
           10 10 10
           11 10 10
           1 1 2
           10 9 10]
n_parts = size(r_parts)[1]
[box[r_parts[i,1],r_parts[i,2],r_parts[i,3]] = mark for i in 1:n_parts]
neigh = [LatticePolymers.n_neighboring_particles(i,r_parts,box,mark,L) for i in 1:n_parts]
@test neigh == [1, 0, 2, 1, 1, 1]

### self_avoiding_cubic_lattice_random_walk ###
n = 10
r = LatticePolymers.self_avoiding_cubic_lattice_random_walk(n)                    
@test size(r) == (n, 3)
# all monomers should have distance 1:
@test [dot(r[i,:]-r[i+1,:],r[i,:]-r[i+1,:]) for i in 1:(n-1)]==ones(Int64,n-1)
# no two monomers should overlap:
for i in 1:(n-1)
    for j in (i+1):n
@test r[i,:]-r[j,:] != [0 0 0]
    end
end 

### self_avoiding_random_walk_in_box ###
nmonos = 20
L = 10

# catch error for L < nmonos
@noinline example() = try
    LatticePolymers.self_avoiding_random_walk_in_box(nmonos, L)
    catch
    return 1
    end
@test example() == 1

L = nmonos+1
box, r = LatticePolymers.self_avoiding_random_walk_in_box(nmonos,L)
@test size(r) == (nmonos, 3)
@test nmonos * 100. - nmonos * 4. - 2. <= sum(box) <= nmonos * 100.

### self_digest_without_attraction ###
L = 3
box=zeros(L,L,L)
N_part=27
m_part=1.
Nsteps = 20
box, r = LatticePolymers.initial_particles_placement(box, L, N_part, m_part)
p_surv = 0. #extremely aggressive self-digest
n_enz = LatticePolymers.self_digest_without_attraction(Nsteps, box, L, r, N_part, m_part, p_surv) 
@test n_enz[end]==1 #nothing should be left except the last enzyme 

### MC_particles_around_polymer_1 ###
nmonos = 10
L = 12
N_parts = 10
m_poly = 1.0
m_part = 200.0
N_steps = 10000
RT = 1.0
E_contact = 0. #no contact energy
box, r = LatticePolymers.self_avoiding_random_walk_in_box(nmonos,L,m_poly,E_contact)
box, r_parts = LatticePolymers.initial_particles_placement(box,L,N_parts,m_part)
energies, particle_contacts = 
    LatticePolymers.MC_particles_around_polymer_1(
        N_steps, box, L, r_parts, N_parts, m_poly, m_part, RT, E_contact
    )
@test mean(particle_contacts) < 0.3

E_contact = -10.0 #strong attraction
box, r = LatticePolymers.self_avoiding_random_walk_in_box(nmonos,L,m_poly,E_contact)
box, r_parts = LatticePolymers.initial_particles_placement(box,L,N_parts,m_part)
energies, particle_contacts = 
    LatticePolymers.MC_particles_around_polymer_1(
        N_steps, box, L, r_parts, N_parts, m_poly, m_part, RT, E_contact
    )
@test mean(particle_contacts) >= 1

# estimate_Boltzmann_weights
N_energies = 100
energies = ones(N_energies)
RT = 1.
E = exp(-1./1.)
Z = N_energies*E
bf = E/Z
@test sum(estimate_Boltzmann_weights(energies, RT) - ones(N_energies).*bf)<1.e-15
