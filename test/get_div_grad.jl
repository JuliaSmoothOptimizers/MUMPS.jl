# Based on Lars Ruthotto's initial implementation.
function get_div_grad(n1::Int, n2::Int, n3::Int)

  # Divergence
  In1 = sparse(1.0I, n1, n1)
  In2 = sparse(1.0I, n2, n2)
  In3 = sparse(1.0I, n3, n3)
  D1 = kron(In3, kron(In2, ddx(n1)))
  D2 = kron(In3, kron(ddx(n2), In1))
  D3 = kron(ddx(n3), kron(In2, In1))

  # DIV from faces to cell-centers
  Div = [D1 D2 D3]

  return Div * Div'
end

# 1D finite difference on staggered grid
ddx(n::Int) = sparse(Bidiagonal(-ones(n), ones(n - 1), :U))
