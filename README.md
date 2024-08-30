# BTTree-Discrete-Free

This package provides tools for deciding whether a subgroup of automorphisms of the Bruhat-Tits tree over a p-adic field is discrete and free, and other algorithms relating to these groups.

This package uses the [BTTree](https://github.com/ariymarkowitz/BTTree) package.

## Main intrinsics

`IsDiscreteFree(T::BTTree, X::[AlgMatElt]: RequireBasis := false) -> BoolElt, .`

Return `false, g` if g is an elliptic element of \<X\>. Otherwise return `true, Y` where Y is a strongly N-reduced basis for \<X\>.

If `RequireBasis := true`, then this function will return `false, g` if g is an element of \<X\> corresponding to a nontrivial word in X that acts trivially on T.

`FundamentalDomain(X::[AlgMatElt], v::BTTVert) -> BTTVert, AlgMatElt`

Return the representative of v in the fundamental domain, with the corresponding group action. X must be a strongly N-reduced basis for \<X\>.

`InFreeIsometryGroup(T::BTTree, g::AlgMatElt, X::[AlgMatElt]) -> BoolElt`

Return true if g is in \<X\>. Assumes X is a strongly N-reduced basis.

`IsSameGroup(T::BTTree, X::[AlgMatElt], Y::[AlgMatElt]) -> BoolElt`

Return true if X and Y generate the same group. Assumes that X and Y are both discrete and free.