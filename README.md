# BTTree-Discrete-Free

This package provides tools for indentifying discrete free subgroups of automorphisms of the Bruhat-Tits tree, and other algorithms relating to these groups.

This package uses the [BTTree](https://github.com/ariymarkowitz/BTTree) package.

## Main intrinsics

`IsDiscreteFree(T::BTTree, X::[AlgMatElt]: RequireBasis := false) -> BoolElt, .`

Return `false, g` if g is an elliptic element of \<X\>. Otherwise return `true, Y` where Y is a reduced free basis of \<X\>.
If `RequireBasis := true`, then this function will return `false, g` if g acts trivially on T.

`FundamentalDomain(X::[AlgMatElt], v::BTTVert) -> BTTVert, AlgMatElt`

Return the representative of v in the fundamental domain, with the corresponding group action.
X must be a strongly N-reduced free basis for \<X\>.

`InFreeIsometryGroup(T::BTTree, g::AlgMatElt, X::[AlgMatElt]) -> BoolElt`

Return true if g is in the group \<X\> with free basis X.

`IsDiscreteFree(T::BTTree, X::[AlgMatElt]: RequireBasis := false) -> BoolElt, .`

Return false, g if g is an elliptic element of \<X\>. Otherwise return true, Y where Y is a reduced free basis of \<X\>.
If RequireBasis := true then this function will return false, g if g acts trivially on T.

`FundamentalDomain(X::[AlgMatElt], v::BTTVert) -> BTTVert, AlgMatElt`

Return the representative of v in the fundamental domain, with the corresponding group action. X must be a strongly N-reduced free basis for \<X\>.

`InFreeIsometryGroup(T::BTTree, g::AlgMatElt, X::[AlgMatElt]) -> BoolElt`

Return true if g is in \<X\>. Assumes X freely generates a discrete free group.

`IsSameGroup(T::BTTree, X::[AlgMatElt], Y::[AlgMatElt]) -> BoolElt`

Return true if g is in the group \<X\> with free basis X. X must be a strongly N-reduced free basis for \<X\>.