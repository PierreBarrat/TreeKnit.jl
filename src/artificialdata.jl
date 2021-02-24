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


function eval_naive_inf(N::Int64, n::Int64, ρ::Float64, simtype::Symbol;
                        Nrep = 1, sfields::Tuple=(:ρ,), out = "")
    #
    args = Dict(:N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype)
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
        mccs = RecombTools.maximal_coherent_clades(t1,t2)
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
    eval_runopt(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol; 
            Tmin=0.001, dT=0.01, Tmax=1.1, M=ceil(n/50), lk_sort=true,
            Nrep = 1,
            sfields::Tuple = ())

Eval the performance of `SplitGraph.runopt` at inferring MCCs. 
"""
function eval_runopt(γ::Real, N::Int64, n::Int64, ρ::Float64, simtype::Symbol; 
    Tmin=0.001, dT=0.01, Tmax=1., Md=10, lk_sort=true,
    Nrep = 1,
    sfields::Tuple = (:ρ,),
    out = "")
    #
    args = Dict(:γ=>γ, :N=>N, :n=>n, :ρ=>ρ, :simtype=>simtype,
        :Md=>Md, :Tmin=>Tmin, :dT=>dT, :Tmax=>Tmax, :lk_sort=>lk_sort)
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
        let γ=γ, Trange=Trange, Md=Md, lk_sort = lk_sort
            t1, t2 = ARGTools.trees_from_ARG(arg)
            mccs = SplitGraph.runopt(t1,t2, γ=γ, Tmin=Tmin, dT=dT, Tmax=Tmax, Md=Md, likelihood_sort=lk_sort)[1]
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

"""
    _eval_mcc_inf(f, N::Int64, n0::Int64, r::Float64; v=true, simtype=:kingman)

Evaluate the inference of MCCs by function `f`. ARG simulation made with `(N,n0,r)`. 
"""
function _eval_mcc_inf(f, N::Int64, n0::Int64, r::Float64; v=true, simtype=:kingman)
    
    arg = ARGTools.SimulateARG.simulate(N, r, n0, simtype=simtype)

    # Real MCCs
    rMCC = ARGTools.MCCs_from_arg(arg); 
    rmMCC = mean([length(x) for x in values(rMCC)])


    # Inferred MCCs
    iMCC = f(arg)
    imMCC = mean([length(x) for x in values(iMCC)])
    t1 = ARGTools.trees_from_ARG(arg)[1]

    # Common branches
    τc_p = 0.; τnc_p = 0.; τc_n = 0.; τnc_n = 0.
    νc_p = 0.; νnc_p = 0.; νc_n = 0.; νnc_n = 0.
    for n in Iterators.filter(x->!x.isroot, values(t1.lnodes))
        if RecombTools.is_branch_in_mccs(n,iMCC) # Branch predicted to be shared with the other tree
            if RecombTools.is_branch_in_mccs(n, rMCC) # Correct prediction
                τc_p += n.data.tau
                νc_p += 1.
            else
                τc_n += n.data.tau
                νc_n += 1.
            end
        else # Branch predicted to not be shared with the other tree
            if !RecombTools.is_branch_in_mccs(n, rMCC) # Correct prediction
                τnc_p += n.data.tau
                νnc_p += 1.
            else
                τnc_n += n.data.tau
                νnc_n += 1.
            end
        end
    end
    νc_p /= (length(t1.lnodes)-1); νnc_p /= (length(t1.lnodes)-1); νc_n /= (length(t1.lnodes)-1); νnc_n /= (length(t1.lnodes)-1)
    T = sum(skipmissing(x.data.tau for x in values(t1.lnodes))) # Total branch length
    τc_p /= T; τnc_p /= T; τc_n /= T; τnc_n /= T

    return Dict(:N=>N, :n0=>n0, :r=>r, 
                :τc_p => τc_p, :τnc_p => τnc_p, :τc_n => τc_n, :τnc_n => τnc_n,
                :νc_p => νc_p, :νnc_p => νnc_p, :νc_n => νc_n, :νnc_n => νnc_n,
                :mMCCsize=>imMCC, :nMCC=>length(iMCC))
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
    if !ismissing(n.data.tau) && n.data.tau < rand(p)
        if !n.isleaf
            nr = delete_node!(n)
            for c in nr.child
                remove_branches!(c, p)
            end
        else
            n.data.tau = 0.
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
"""
    remove_branches!(t::Tree, p) 
"""
function remove_branches(t::Tree, p)
    tt = deepcopy(t)
    remove_branches!(tt, p)
    node2tree!(tt, tt.root)
    return tt
end