intrinsic Length(T::BTTree, g::AlgMatElt) -> RngIntElt
{ Return the length of g }
  v := Origin(T);
  return Distance(v, v*g);
end intrinsic;

intrinsic LexicalCmp(T::BTTree, v::BTTVert, w::BTTVert) -> RngIntElt
{ Compare the paths from the origin to v and w respectively under the lexical order }
  if v eq w then return 0; end if;
  o := Origin(T);
  p := Path(o, v);
  q := Path(o, w);

  return LexicalCmp(p, q);
end intrinsic;

intrinsic LexicalCmp(p::[BTTVert], q::[BTTVert]) -> RngIntElt
{ Return -1, 1 or 0 if p < q, p > q, or p = q respectively under the lexical order }
  if (#p lt #q) then
    return -1;
  elif (#q lt #p) then
    return 1;
  end if;

  i := 1;
  for i in [1..#p] do
    if p[i] ne q[i] then
      pEdge := TypeOfNeighbor(p[i-1], p[i]);
      qEdge := TypeOfNeighbor(q[i-1], q[i]);
      return pEdge lt qEdge select -1 else 1;
    end if;
  end for;

  return 0;
end intrinsic;

intrinsic LexicalCmp(T::BTTree, g::AlgMatElt, h::AlgMatElt) -> RngIntElt
{ Return -1, 1 or 0 if g < h, g > h, or g ≃ h respectively under the lexical preorder }
  v := Origin(T);
  p := Path(v, v*g);
  q := Path(v, v*h);

  if (#p lt #q) then
    return -1;
  elif (#q lt #p) then
    return 1;
  end if;

  pInv := Path(v, v*g^-1);
  qInv := Path(v, v*h^-1);

  // Take the initial half of each path (Note: We halve the number of edges, not vertices).
  len := (#p+1) div 2;
  pHalf := p[1..len];
  qHalf := q[1..len];
  pInvHalf := pInv[1..len];
  qInvHalf := qInv[1..len];

  if LexicalCmp(pHalf, pInvHalf) eq -1 then
    pMin := pHalf;
    pMax := pInvHalf;
  else
    pMin := pInvHalf;
    pMax := pHalf;
  end if;

  if LexicalCmp(qHalf, qInvHalf) eq -1 then
    qMin := qHalf;
    qMax := qInvHalf;
  else
    qMin := qInvHalf;
    qMax := qHalf;
  end if;

  cmp := LexicalCmp(pMin, qMin);
  if cmp eq 0 then return LexicalCmp(pMax, qMax); end if;
  return cmp;
end intrinsic;

intrinsic ReduceGenerators(T::BTTree, X::[AlgMatElt]) -> RngIntElt, .
{ Return 0, i if i is the index of an elliptic in X.
If X is not reduced, return 1, Y where Y is a reduction of X via a Nielsen transformation.
Otherwise return 3. }
  // Check for elliptic isometries.
  for i -> g in X do
    if IsElliptic(g) then
      return 0, i;
    end if;
  end for;

  // Check for a smaller isometry.
  for i -> g in X do
    for j -> h in X do
      if i eq j then continue; end if;

      // Since X is a sequence, this automatically handles a violation of A1.
      for str in ["h*g", "g*h^-1", "g*h", "h^-1*g"] do
        a := eval(str);
        if LexicalCmp(T, g, a) eq 1 then
          return 1, Insert(X, i, i, [a]);
        end if;
      end for;
    end for;
  end for;

  return 3, _;
end intrinsic;

intrinsic IsDiscreteFree(T::BTTree, X::[AlgMatElt]: RequireBasis := false) -> BoolElt, .
{ Return false, g if g is an elliptic element of <X>. Otherwise return true, Y where Y is a reduced free basis of <X>.
If RequireBasis := true then this function will return false, g if g acts trivially on T. }
  Y := X;
  while true do
    a, b := ReduceGenerators(T, Y);
    if a eq 0 then
      // b is the index of a scalar matrix.
      if not RequireBasis then
        // Check that Y approximates a scalar matrix.
        A := Y[b];
        if IsApproximatelyIdentity(A) then
          Remove(~Y, b);
          continue;
        end if;
      end if;
      return false, Y[b];
    elif a eq 1 then
      Y := b;
      continue;
    end if;
    return true, Y;
  end while;
end intrinsic;

intrinsic FundamentalDomain(v::BTTVert, X::[AlgMatElt]) -> BTTVert, AlgMatElt
{ Return the representative of v in the fundamental domain, with the corresponding group action.
 X must be a strongly N-reduced free basis for <X>. }
  T := Parent(v);
  g := MatrixAlgebra(Field(v), 2)!1;
  w := v;
  done := true;
  while done eq true do
    done := false;
    for h in X do
      if LexicalCmp(T, w, w*h) eq 1 then
        w := w * h;
        g := g * h;
        done := true;
        break;
      elif LexicalCmp(T, w, w*h^-1) eq 1 then
        w := w * h^-1;
        g := g * h^-1;
        done := true;
        break;
      end if;
    end for;
  end while;
  return w, g;
end intrinsic;

intrinsic InFreeIsometryGroup(T::BTTree, g::AlgMatElt, X::[AlgMatElt]) -> BoolElt
{ Return true if g is in the group <X> with free basis X }
  _, h := FundamentalDomain(Origin(T)*g^(-1), X);
  return h eq g;
end intrinsic;

intrinsic IsSameGroup(T::BTTree, X::[AlgMatElt], Y::[AlgMatElt]) -> BoolElt
{ Return true if X and Y generate the same group. Assumes that X and Y are both discrete and free. }
  b1, X2 := IsDiscreteFree(T, X);
  b2, Y2 := IsDiscreteFree(T, Y);
  if not (b1 and b2) then return false; end if;

  cmp := func<g, h | LexicalCmp(T, g, h)>;
  X2 := Sort(X2, cmp);
  Y2 := Sort(Y2, cmp);

  if #X2 ne #Y2 then
    return false;
  end if;

  for i -> x in X2 do
    if not (IsApproximatelyIdentity(x*Y2[i]) or IsApproximatelyIdentity(x*Y2[i]^-1)) then
      return false;
    end if;
  end for;

  return true;
end intrinsic