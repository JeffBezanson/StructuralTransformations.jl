module __StructuralTransformations

function __init__()
    restored = Base._require_from_serialized(joinpath(@__DIR__, "qhONu_9ZVWO.ji"))
    if restored isa Exception
        throw(restored)
    end
    global StructuralTransformations = restored[end]
end

end
