% A wrapper function that makes sure the cooling function fits the problem
% description.
function [func]=getBoundaryFunc(f)
    assert(f(0)==980)
    func=@(t) max(f(t), 20);
end