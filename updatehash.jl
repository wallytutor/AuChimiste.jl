# -*- coding: utf-8 -*-
using Pkg.GitTools
using SHA
using TOML

TREE_HASH = bytes2hex(Pkg.GitTools.tree_hash(joinpath(@__DIR__, "artifacts")))

function file_hash(artifacts, fpath)
    base = basename(fpath)
    name = join(split(base, ".")[1:end-1], "_")
    hash = bytes2hex(sha256(read(fpath)))
    
    artifacts[name] = Dict(
        "files"         => [base], 
        "hash"          => hash,
        "git-tree-sha1" => TREE_HASH
    )
    
    @info("""\
    Generated SHA-256 hash for $(fpath)
    
    Using SHA $(hash) for artifact name `$(name)`
    """)
end

function update_sha256(; artifacts, files)
    artifacts_data = TOML.parsefile(artifacts)

    for fpath in files
        file_hash(artifacts_data, joinpath(@__DIR__, "artifacts", fpath))
    end

    open(artifacts, "w") do fp
        TOML.print(fp, artifacts_data)
    end
end

update_sha256(; artifacts = joinpath(@__DIR__, "Artifacts.toml"), 
                files     = [
                    "thermo_compounds.yaml",
                ])
