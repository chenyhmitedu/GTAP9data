module GTAP9data

using CSV, DataFrames, JLD2

# Read data into a Vector
path        = "./src/Input_GTAP9/"
csv_files   = filter(f -> endswith(f, ".csv"), readdir(path))
n           = length(csv_files)
dic_df      = Dict()
dic_tmp     = Dict()
d     = Dict()
s     = Dict()

# Convert dic_df (key => dataframe) into d (key => dictionary) for a dataframe with # column >= 2, and into s (vector) with # column < 2
for i ∈ 1:n
    str_tmp = csv_files[i][1:end-4]
    dic_df[str_tmp]    = CSV.read(joinpath(path, csv_files[i]), DataFrame, delim=",")
    cols    = names(dic_df[str_tmp])
    
    if ncol(dic_df[str_tmp]) > 2
        dic_tmp = Dict(
                Symbol.(Tuple(row[cols[1:end-1]])) => row[cols[end]]
                for row in eachrow(dic_df[str_tmp])
                )
        d[str_tmp]    = dic_tmp

    elseif ncol(dic_df[str_tmp]) == 2
        dic_tmp = Dict(
                Symbol.(row[cols[1]]) => row[cols[2]]
                for row in eachrow(dic_df[str_tmp])
                )
        d[str_tmp]    = dic_tmp

    else
        s[str_tmp]    = Symbol.(Vector(dic_df[str_tmp][:, 1]))
    end
end

# Fill up each missing key element with a (k, 0) 
function fullspace(dict::Dict, args...)
    key_space = collect(Iterators.product(args...))
    for k in key_space
        if !haskey(dict, k)
            dict[k] = 0
        end
    end
    return dict
end

function fullspace(dict::Dict, args)
    key_space = collect(args)
    for k in key_space
        if !haskey(dict, k)
            dict[k] = 0
        end
    end
    return dict
end

# Expand to the full set space and use GTAP notation
vdfm    = fullspace(d["vdfm"], s["set_i"], s["set_g"], s["set_r"])
vxmd    = fullspace(d["vxmd"], s["set_i"], s["set_r"], s["set_r"])
vst     = fullspace(d["vst"], s["set_i"], s["set_r"])
rtms    = fullspace(d["rtms"], s["set_i"], s["set_r"], s["set_r"])
rtxs    = fullspace(d["rtxs"], s["set_i"], s["set_r"], s["set_r"])
vifm    = fullspace(d["vifm"], s["set_i"], s["set_g"], s["set_r"])
rtfd    = fullspace(d["rtfd"], s["set_i"], s["set_g"], s["set_r"])
rtfi    = fullspace(d["rtfi"], s["set_i"], s["set_g"], s["set_r"])
rto     = fullspace(d["rto"], s["set_g"], s["set_r"])
vfm     = fullspace(d["vfm"], s["set_f"], s["set_g"], s["set_r"])
rtf     = fullspace(d["rtf"], s["set_f"], s["set_g"], s["set_r"])
vtwr    = fullspace(d["vtwr"], s["set_i"], s["set_i"], s["set_r"], s["set_r"])
esubd   = fullspace(d["esubd"], s["set_i"])
esubm   = fullspace(d["esubm"], s["set_i"])
esubva  = fullspace(d["esubva"], s["set_g"])
set_i   = s["set_i"]
set_g   = s["set_g"]
set_r   = s["set_r"]
set_f   = s["set_f"]
set_cgi = setdiff(set_g, set_i)

# Assignment done in GTAPinGAMS
d["esub"]       = Dict(i => 0 for i ∈ s["set_g"])       # Top-level elasticity of substitution
d["esub"][:c]   = 1                                     

vdm = Dict(
    (i, r) => sum(vdfm[(i, g, r)] for g in set_g)
    for i in set_i, r in set_r
)

vom = Dict(
    (i, r) => sum(vxmd[(i, r, s)] for s in set_r) + vdm[(i, r)] + vst[(i, r)]
    for i ∈ set_i, r ∈ set_r
)

vom_a = Dict(
    (i, r) => sum((vdfm[(j, i, r)]*(1+rtfd[(j, i, r)]) + vifm[(j, i, r)]*(1+rtfi[(j, i, r)]))/(1-rto[(i, r)]) for j ∈ set_i) 
    for i ∈ set_cgi, r ∈ set_r
)

vom = merge(vom, vom_a)

vdm_a = Dict(
    (i, r) => vom[(i, r)] for i ∈ [:c, :g], r ∈ set_r
)

vdm = merge(vdm, vdm_a)

#test = sum(vom[(i, r)] for i ∈ set_g, r ∈ set_r)
#key = [k for (k, v) in vom if v == 0]
#key = [k for (k, v) in vdm if v == 0]
#key = [k for (k, v) in vtw if v == 0]
#key = [k for (k, v) in vim if v == 0]

pvxmd = Dict(
    (i, s, r) => (1 + rtms[(i, s, r)])*(1 - rtxs[(i, s, r)]) for i ∈ set_i, s ∈ set_r, r ∈ set_r
)

pvtwr = Dict(
    (i, s, r) => 1 + rtms[(i, s, r)] for i ∈ set_i, s ∈ set_r, r ∈ set_r
)

vtw = Dict(
    j => sum(vst[(j, r)] for r ∈ set_r)
    for j ∈ set_i
)

vim = Dict(
    (i, r) => sum(vifm[(i, g, r)] for g ∈ set_g)
    for i ∈ set_i, r ∈ set_r
)

vb = Dict(
    r => sum(vom[(i, r)] for i ∈ set_cgi)
        -sum(vfm[(f, g, r)] for f ∈ set_f, g ∈ set_g)
        -sum(vom[(g, r)]*rto[(g, r)] for g ∈ set_g)
        -sum(vdfm[(i, g, r)]*rtfd[(i, g, r)] + vifm[(i, g, r)]*rtfi[(i, g, r)] for i ∈ set_i, g ∈ set_g)
        -sum(vfm[(f, g, r)]*rtf[(f, g, r)] for f ∈ set_f, g ∈ set_g)
        -sum(rtms[(i, s, r)]*(vxmd[(i, s, r)]*(1-rtxs[(i, s, r)]) + sum(vtwr[(j, i, s, r)] for j ∈ set_i)) for i ∈ set_i, s ∈ set_r)
        +sum(rtxs[(i, r, s)]*vxmd[(i, r, s)] for i ∈ set_i, s ∈ set_r)
    for r ∈ set_r
)

end # module GTAP9data
