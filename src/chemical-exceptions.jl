# -*- coding: utf-8 -*-

export NoSuchElementError
export NoIsotopeProvidedError
export EmptyCompositionError
export InvalidScalerError

"Element (or isotope) was not found in user database."
struct NoSuchElementError <: ChemicalElementsError
    message::String

    function NoSuchElementError(e)
        return new("""\
            No such element $(e) in the elements dictionary. If you are \
            trying to access an isotope, please make sure you create it \
            before.
            """)
    end
end

"Unstable elements do not provide atomic mass."
struct NoIsotopeProvidedError <: ChemicalElementsError
    message::String

    function NoIsotopeProvidedError(e)
        return new("""\
            Accessing the atomic mass of unstable element $(e) is not \
            supported. Please consider creating a named isotope of \
            this element with `add_isotope`.
            """)
    end
end

"A composition set is missing for the given component."
struct EmptyCompositionError <: ChemicalComponentsError
    message::String

    function EmptyCompositionError()
        return new("""\
            Cannot create a composition with an empty set. Please \
            supply the elements and related amounts accoding to the \
            type of composition specification.
            """)
    end
end

"The provided scaler targets an unspecified element."
struct InvalidScalerError <: ChemicalComponentsError
    message::String

    function InvalidScalerError(s, e)
        return new("""\
            Scaling configuration must be a valid element from the \
            system. Could not find $(s) among $(e).
            """)
    end
end
