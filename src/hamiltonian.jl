struct AbelianSchwingerModel
    nsites::Int
    nflavors::Int
    m::Float64
    g::Float64
    a::Float64
    l0::Int  # background field
    cutoff::Float64
    ham::Vector{MPO}  # dimensionless
end

function AbelianSchwingerModel(nsites::Int, nflavors::Int, m::Float64, g::Float64, a::Float64, l0::Int=0; cutoff=1e-9)
    @assert nsites > 0
    @assert nflavors > 0
    @assert m > 0
    @assert a > 0
    ham = if nflavors == 1
        MPO_AbelianSchwingerSF(nsites, g, m, a, l0; cutoff)
    else
        MPO_AbelianSchwingerMF(nsites, nflavors, g, m, a, l0; cutoff)
    end

    return AbelianSchwingerModel(nsites, nflavors, m, g, a, l0, cutoff, ham)
end


"""Get the MPO for the single flavor U(1) Schwinger model Hamiltonian"""
function MPO_AbelianSchwingerSF(
            nsites::Int,
            m::Float64,
            g::Float64,
            a::Float64=1.0,
            l0::Int=0;
            cutoff::Float64=1e-9)
    
    os = OpSum()
    # dimensionless parameters
    x = inv(a^2 * g^2)
    μ = 2m * inv(a * g^2)

    sites = siteinds("S=1/2", nsites)
    for i in 1:nsites-1
        # hopping term
        os += x, "S+", i, "S-", i + 1
        os += x, "S-", i, "S+", i + 1

        # mass term
        os += 0.5μ, "Id", i
        os += 0.5μ * (-1)^(i - 1), "Sz", i
    end

    # mass terms for last site
    os += 0.5μ, "Id", nsites
    os += 0.5μ * (-1)^(i - 1), "Sz", nsites
    W_SF = MPO(os, sites)

    # Gauss law encoding term
    for i in 1:nsites-1
        os_gi = OpSum()
        os_gi += l0, "Id", 1  # background field
        for j in 1:i
            os_gi += 0.5 * (-1)^(k - 1), "Id", j
            os_gi += 0.5, "Sz", j
        end
        W_gi = MPO(os_hi, sites)
        W_gi = apply(W_gi, W_gi; cutoff)

        W_SF = add(W_0, W_gi; cutoff)
    end

    # enforce zero-net charge i.e half-filling through a strong penalty on net spin
    os_p = OpSum()
    for i in 1:nsites
        os += 1000μ, "Sz", i
    end
    W_penalty = MPO(os_p, sites)
    W_penalty = apply(W_penalty, W_penalty; cutoff)  # (Sum of Sz)^2 ~ 0

    return [W_SF, W_penalty]
end

function MPO_AbelianSchwingerMF(
            nsites::Int,
            nflavors::Int,
            m::Float64,
            g::Float64,
            a::Float64=1.0,
            l0::Int=0;
            cutoff::Float64=1e-9)
    
    return nothing
end