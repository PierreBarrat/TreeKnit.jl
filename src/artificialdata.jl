export eval_mcc_inf, df_fields, eval_runopt
export remove_branches!

#####################
##################### Evaluating inference of splits
#####################




#####################
##################### Evaluating MCC inference
#####################

function get_r(ρ::Float64, n::Int64, N::Int64, simtype::Symbol)
    if simtype == :kingman
        return ρ * n / N
    elseif simtype == :yule
        return ρ / N
    else
        @error "Unrecognized `simtype`."
    end
end

"""
    eval_naive_inf(N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
        cutoff = 0.,
        Nrep = 1,
        sfields::Tuple = (:ρ,:cutoff),
        out = ""
    )
"""
function eval_naive_inf(N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
    cutoff = 0.,
    Nrep = 1,
    sfields::Tuple = (:ρ,:cutoff),
    out = ""
)
    #
    args = Dict(:N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype, :cutoff=>cutoff)
    #
    dat = DataFrame(df_fields())
    for f in sfields
        insertcols!(dat, f=>Any[])
    end
    #
    r = get_r(ρ, n, N, simtype)
    #
    function f(arg::ARGTools.ARG)
        t1, t2 = ARGTools.trees_from_ARG(arg)
        if cutoff != 0
            remove_branches!(t1, Distributions.Exponential(cutoff*N))
            remove_branches!(t2, Distributions.Exponential(cutoff*N))
        end
        init_splits = Dict(1=>SplitList(t1), 2=>SplitList(t2))
        trees = Dict(1=>t1, 2=>t2)
        MCCs, resolved_splits = computeMCCs!(trees, naive=true)
        return MCCs, resolved_splits, init_splits
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r, simtype=simtype)
        d = td[1]
        args[:time] = td[2]
        args[:bytes] = td[3]
        for f in sfields
            d[f] = args[f]
        end
        push!(dat, d)
    end
    #
    if out != ""
        CSV.write(out, dat)
    end
    #
    return dat
end

"""
    eval_real(N::Int, n::Int, ρ::Float64, simtype::Symbol;
        Nrep = 1,
        out = "",
    )
"""
function eval_real(N::Int, n::Int, ρ::Float64, simtype::Symbol;
    Nrep = 1,
    out = "",
)
    dat = DataFrame(df_fields())
    insertcols!(dat, :ρ=>Float64[])
    r = get_r(ρ, n, N, simtype)
    #
    function f(arg::ARGTools.ARG)
        t1, t2 = ARGTools.trees_from_ARG(arg)
        init_splits = Dict(1=>SplitList(t1), 2=>SplitList(t2))
        resolved_splits = Dict(1=>[], 2=>[])
        MCCs = Dict((1,2) => ARGTools.MCCs_from_arg(arg))
        return MCCs, resolved_splits, init_splits
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r, simtype=simtype)
        d = td[1]
        d[:ρ] = ρ
        push!(dat, d)
    end
    #
    if out != ""
        CSV.write(out, dat)
    end
    #
    return dat
end


"""
    eval_runopt(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
    Tmin=0.001, dT=0.01, Tmax=1., Md=10, lk_sort=true,
    cutoff = 0.,
    Nrep = 1,
    preresolve = true,
    crossmap_prune=false,
    crossmap_resolve=false,
    simulate_seqs = false,
    sfields::Tuple = (:ρ,:cutoff, :preresolve, :crossmap_prune, :crossmap_resolve),
    out = "",
    verbose=false
)

Eval the performance of `SplitGraph.runopt` at inferring MCCs.
"""
function eval_runopt(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
    Tmin=0.001, dT=0.01, Tmax=1., Md=10, lk_sort=true,
    cutoff = 0.,
    Nrep = 1,
    preresolve = true,
    crossmap_resolve = false,
    crossmap_prune = false,
    suspmut_threshold = 2,
    simulate_seqs = false,
    sfields::Tuple = (:ρ,:cutoff, :preresolve, :crossmap_prune, :crossmap_resolve),
    out = "",
    verbose=false
)
    #
    args = Dict(:γ=>γ, :N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype,
        :Md=>Md, :Tmin=>Tmin, :dT=>dT, :Tmax=>Tmax, :lk_sort=>lk_sort,
        :cutoff=>cutoff, :preresolve=>preresolve,
        :crossmap_resolve=>crossmap_resolve, :crossmap_prune=>crossmap_prune,
        :suspmut_threshold=>suspmut_threshold)
    #
    dat = DataFrame(df_fields())
    for f in sfields
        insertcols!(dat, f=>Any[])
    end
    #
    r = get_r(ρ, n, N, simtype)
    Trange = reverse(Tmin:dT:Tmax)
    #
    function f(arg::ARGTools.ARG)
        let cutoff=cutoff, N=N
            t1, t2 = ARGTools.trees_from_ARG(arg)
            trees = Dict(1=>t1, 2=>t2)
            if cutoff != 0 && !crossmap_prune && !crossmap_resolve
                remove_branches!(t1, Distributions.Exponential(cutoff*N))
                remove_branches!(t2, Distributions.Exponential(cutoff*N))
            elseif crossmap_prune || crossmap_resolve || simulate_seqs
                simulate_sequences!(trees, N, cutoff)
            end
            init_splits = Dict(1=>SplitList(t1), 2=>SplitList(t2))
            oa = OptArgs(
                γ=γ,
                Tmin=Tmin, dT=dT, Tmax=Tmax,
                Md=Md,
                likelihood_sort=lk_sort,
                crossmap_resolve=crossmap_resolve,
                crossmap_prune=crossmap_prune,
                suspmut_threshold=suspmut_threshold,
                verbose=verbose,
            )
            try
                MCCs, resolved_splits = computeMCCs!(trees, oa, preresolve=preresolve)
                return MCCs, resolved_splits, init_splits
            catch err
                mkpath("tmp")
                write_newick("tmp/tree1.nwk", trees[1])
                write_newick("tmp/tree2.nwk", trees[2])
                #TreeTools.write_fasta("tmp/aln1.fasta", trees[1], :selfseq)
                #TreeTools.write_fasta("tmp/aln2.fasta", trees[2], :selfseq)
                error(err)
            end
        end
    end
    #
    for rep in 1:Nrep
        td = @timed _eval_mcc_inf(f, N, n, r, simtype=simtype)
        d = td[1]
        args[:time] = td[2]
        args[:bytes] = td[3]
        for f in sfields
            d[f] = args[f]
        end
        push!(dat, d)
    end
    #
    if out != ""
        CSV.write(out, dat)
    end
    #
    return dat
end


function simulate_sequences!(trees::Dict, N, c = 0.5, L = 1000; polytomies=true)
    μ = c == 0 ? 1 / (N*L*0.2) * 1.25 : 1 / (N*L*c) * 1.25
    # Simulate sequences && prune branches with no mutations
    for (i,t) in trees
        TreeAlgs.Evolve.evolve!(t, L, μ;  seqkey = :selfseq)
        TreeTools.compute_mutations!(n -> n.data.dat[:selfseq], t, :realmuts)
        c != 0 && polytomies && TreeTools.delete_branches!(n->isempty(n.data.dat[:realmuts]), t)
    end
    # Cross-map sequences
    crossmap_sequences!(trees, :selfseq, :cmseq)
    #
    return nothing
end

function _eval_mcc_inf(rMCC, iMCC, t1::Tree)
    # Mean MCC sizes
    rmMCC = mean([length(x) for x in values(rMCC)])
    imMCC = mean([length(x) for x in values(iMCC)])
    # Common branches
    τc_p = 0.; τnc_p = 0.; τc_n = 0.; τnc_n = 0.
    νc_p = 0.; νnc_p = 0.; νc_n = 0.; νnc_n = 0.
    for n in Iterators.filter(x->!x.isroot, values(t1.lnodes))
        if RecombTools.is_branch_in_mccs(n,iMCC) # Branch predicted to be shared with the other tree
            if RecombTools.is_branch_in_mccs(n, rMCC) # Correct prediction
                τc_p += n.tau
                νc_p += 1.
            else
                τc_n += n.tau
                νc_n += 1.
            end
        else # Branch predicted to not be shared with the other tree
            if !RecombTools.is_branch_in_mccs(n, rMCC) # Correct prediction
                τnc_p += n.tau
                νnc_p += 1.
            else
                τnc_n += n.tau
                νnc_n += 1.
            end
        end
    end
    νc_p /= (length(t1.lnodes)-1); νnc_p /= (length(t1.lnodes)-1); νc_n /= (length(t1.lnodes)-1); νnc_n /= (length(t1.lnodes)-1)
    T = sum(skipmissing(x.tau for x in values(t1.lnodes))) # Total branch length
    τc_p /= T; τnc_p /= T; τc_n /= T; τnc_n /= T

    return Dict(:τc_p => τc_p, :τnc_p => τnc_p, :τc_n => τc_n, :τnc_n => τnc_n,
                :νc_p => νc_p, :νnc_p => νnc_p, :νc_n => νc_n, :νnc_n => νnc_n,
                :mMCCsize=>imMCC, :nMCC=>length(iMCC))
end

function _eval_split_inf(true_splits, resolved_splits, init_splits)
    rs_p = 0
    rs_n = 0
    rs_init = length(init_splits) / length(true_splits)
    for rs in resolved_splits
        if in(rs, true_splits, usemask=false)
            rs_p += 1
        else
            rs_n += 1
        end
    end
    return Dict(:rs_p => rs_p / length(true_splits), :rs_n=> rs_n / length(true_splits),
                :rs_init => rs_init)

end

"""
    _eval_mcc_splits_inf(f, N::Int64, n0::Int64, r::Float64; v=true, simtype=:kingman)

Evaluate the inference of MCCs by function `f`. ARG simulation made with `(N,n0,r)`.
"""
function _eval_mcc_inf(f, N::Int64, n0::Int64, r::Float64; v=true, simtype=:kingman)

    arg = ARGTools.SimulateARG.simulate(N, r, n0, simtype=simtype)

    # Real MCCs
    rMCC = ARGTools.MCCs_from_arg(arg);
    trees = ARGTools.trees_from_ARG(arg)
    true_splits = Dict(i=>SplitList(t) for (i,t) in enumerate(trees))

    # Inferred MCCs
    iMCC, resolved_splits, init_splits = f(arg)
    t1 = trees[1]

    out = _eval_mcc_inf(rMCC, iMCC[1,2], t1)
    out_ = _eval_split_inf(true_splits[1], resolved_splits[1], init_splits[1])
    for (k,v) in out_
        out[k] = v
    end
    out[:N] = N
    out[:n0] = n0
    out[:r] = r
    return out
end



function eval_mcc_inf(f, N::Int64, n0::Int64, rrange;
        Nrep=25*ones(Int64, length(rrange)), v=true, simtype=:kingman)
    dat = DataFrame(df_fields())
    for (i,r) in enumerate(rrange)
        v && print("r=$(round(ρ,digits=5)), $i/$(length(rrange))                                                     \r")
        Nr = Nrep[i]
        for rep in 1:Nr
            push!(dat, _eval_mcc_inf(f, N, n0, r, v=v, simtype=simtype))
        end
    end
    return dat
end
function eval_mcc_inf(f, N::Int64, n0::Int64, r::Float64;
        Nrep=25, v=true, simtype=:kingman)
    dat = DataFrame(df_fields())
    for rep in 1:Nrep
        push!(dat, _eval_mcc_inf(f, N, n0, r, v=v, simtype=simtype))
    end
    return dat
end

"""
    df_fields()

Fields returned by `_eval_mcc_inf`:
- `N`, `n0`, `r`
- `ν(n)c_p(n)`: Fraction of inferred (non-)common branches, true (p) or false (n)
- `τ(n)c_p(n)`: Length of inferred (non-)common branches, true (p) or false (n) (normalized)
- `mMCCsize`: mean inferred MCC size
- `nMCC`: number of inferred MCCs
"""
function df_fields()
    return (N = Int64[], n0 = Int64[], r = Float64[], # Simulation parameters
        νc_p = Float64[], νnc_p = Float64[], # Fraction of correctly inferred common and non-common branches
        νc_n = Float64[], νnc_n = Float64[], # Fraction of incorrectly inferred common and non-common branches
        τc_p = Float64[], τnc_p = Float64[], # Total length of correctly inferred common and non-common branches
        τc_n = Float64[], τnc_n = Float64[], # Total length of incorrectly inferred common and non-common branches
        rs_p = Float64[], rs_n = Float64[], # Number of (in)-correctly resolved splits, scaled by total number of true splits
        rs_init = Float64[], # Initial number of splits scaled by total number of true splits
        mMCCsize = Float64[], # Mean MCC size (inferred)
        nMCC = Int64[], # Number of MCCs (inferred)
        )
end


#####################
##################### REMOVING BRANCHES
#####################
"""
    remove_branches!(n::TreeNode, p::Distribution)
"""
function remove_branches!(n::TreeNode, p::Distribution)
    if !ismissing(n.tau) && n.tau < rand(p)
        if !n.isleaf
            nr = delete_node!(n)
            for c in nr.child
                remove_branches!(c, p)
            end
        else
            n.tau = 0.
        end
    else
        for c in n.child
            remove_branches!(c, p)
        end
    end
end

"""
    remove_branches!(t::Tree, p)

Stochastically remove branches from `t` using probability distribution `p`
"""
function remove_branches!(t::Tree, p)
    remove_branches!(t.root, p)
    node2tree!(t, t.root)
    nothing
end
remove_branches!(t::Tree, c::Real, N::Int) = remove_branches!(t, Exponential(c*N))


