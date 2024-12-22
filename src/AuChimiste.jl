# -*- coding: utf-8 -*-
module AuChimiste

function __init__()
    # Path to annex modules (add it only if not already there)
    let CHIMISTE_PATH = abspath(@__DIR__)
        (CHIMISTE_PATH in LOAD_PATH) && return nothing
        push!(LOAD_PATH, CHIMISTE_PATH)
    end

    return nothing
end

end # (module AuChimiste)