function [func]=getBoundaryFunc(f)
    assert(f(0)==980)
    func=@(t) max(f(t), 20);
end